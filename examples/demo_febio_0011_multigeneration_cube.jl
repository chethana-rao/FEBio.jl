using Comodo # https://github.com/COMODO-research/Comodo.jl
using FEBio
using Printf
using GLMakie
using LinearAlgebra
using XML
using GeometryBasics

######
# Set FEBio path here for now
const FEBIO_PATH = "/home/simulimb/FEBioStudio/bin/febio4" # Path to FEBio executable

###### 
# Control parameters 
sampleSize = 10.0
pointSpacing = 1.0

#Material Constants
c1 = 2 #Shear-modulus-like parameter
m1 = 2  #Material parameter setting degree of non-linearity
k_factor = 50  #Bulk modulus factor 
k = c1*k_factor #Bulk modulus
c1_g = [c1/1000 c1*2]
k_g = c1_g.*k_factor
###### 
# Creating a hexahedral mesh for a cube 
boxDim = sampleSize.*[1,1,1] # Dimensionsions for the box in each direction
boxEl = ceil.(Int64,boxDim./pointSpacing) # Number of elements to use in each direction 
E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)

# Create face sets to define node sets later 
Fb_top = Fb[CFb_type.==2]
Fb_bottom = Fb[CFb_type.==1]
Fb_s1 = Fb[CFb_type.==6]
Fb_s2 = Fb[CFb_type.==3]
indNodes_bottom = elements2indices(Fb_bottom)
indNodes_top = elements2indices(Fb_top)
######
# Defining displacement of the top surface in terms of x, y, and z components
strainApplied = 0.3 # Equivalent linear strain
loadingOption ="tension" # "tension" or "compression"

if loadingOption=="tension"
    stretchLoad = 1 + strainApplied
elseif loadingOption=="compression"
    stretchLoad = 1- strainApplied
end

displacementMagnitude = (stretchLoad*sampleSize) - sampleSize 
#Splitting the material into two domains 
VE = simplexcenter(E,V)
XE = [v[1] for v in VE]
X = [v[1] for v in V]
Xmax = maximum(X)
Xmin = minimum(X)
logicMaterial_1 = XE .<= (Xmax+Xmin)/2
logicMaterial_2 = logicMaterial_1 .==0
indshow_1 = findall(logicMaterial_1)
indshow_2 = findall(logicMaterial_2)
Fs = element2faces(E[indshow_2])
Fs,Vs = separate_vertices(Fs,V)
Ms = GeometryBasics.Mesh(Vs,Fs)

Em1= E[indshow_1]
Em2= E[indshow_2]
# Define file names
saveDir = joinpath(febiojl_dir(),"assets","temp") # Main directory to save FEBio input and output files
if !isdir(saveDir)
    mkdir(saveDir)      
end

filename_FEB = joinpath(saveDir,"febioInputFile_01.feb")   # The main FEBio input file
filename_xplt = joinpath(saveDir,"febioInputFile_01.xplt") # The XPLT file for viewing results in FEBioStudio
filename_log = joinpath(saveDir,"febioInputFile_01_LOG.txt") # The log file featuring the full FEBio terminal output stream
filename_disp = "febioInputFile_01_DISP.txt" # A log file for results saved in same directory as .feb file  e.g. nodal displacements
filename_energy_out = "febioInputFile_energy_out.txt"
filename_stress = "febioInputFile_01_STRESS.txt"
######
# Define febio input file XML
doc,febio_spec_node = feb_doc_initialize()

aen(febio_spec_node,"Module"; type = "solid") # Define Module node: <Module type="solid"/>


Globals_node   = aen(febio_spec_node,"Globals")

Constants_node = aen(Globals_node,"Constants")
    aen(Constants_node,"R",8.3140000e-06)
    aen(Constants_node,"T",298)
    aen(Constants_node,"F",9.6485000e-05)

Material_node = aen(febio_spec_node,"Material")

material_node = aen(Material_node,"material"; id = 1, name="Normal_material", type="Ogden unconstrained")
    aen(material_node,"c1",c1)
    aen(material_node,"m1",m1)
    aen(material_node,"c2",c1)
    aen(material_node,"m2",-m1)
    aen(material_node,"cp",k)

    material_node = aen(Material_node,"material"; id = 2, name="Multigen_material", type="multigeneration")
    generation_node = aen(material_node, "generation";id = 1)
    aen(generation_node, "start_time",0)
    solid_node = aen(generation_node,"solid"; type = "Ogden unconstrained")
    aen(solid_node,"c1",c1_g[1])
    aen(solid_node,"m1",m1)
    aen(solid_node,"c2",c1_g[1])
    aen(solid_node,"m2",-m1)
    aen(solid_node,"cp",k_g[1])
    generation_node = aen(material_node, "generation";id = 2)
    aen(generation_node, "start_time",1)
    solid_node = aen(generation_node,"solid"; type = "Ogden unconstrained")
    aen(solid_node,"c1",c1_g[2])
    aen(solid_node,"m1",m1)
    aen(solid_node,"c2",c1_g[2])
    aen(solid_node,"m2",-m1)
    aen(solid_node,"cp",k_g[2])


Mesh_node = aen(febio_spec_node,"Mesh")

Nodes_node = aen(Mesh_node,"Nodes"; name="nodeSet_all")
    for q ∈ eachindex(V)
        # aen(Nodes_node,"node",@sprintf("%.2f, %.2f, %.2f",V[q][1],V[q][2],V[q][3]); id = q)
        aen(Nodes_node,"node", join([@sprintf("%.16e",x) for x ∈ V[q]],','); id = q)     
    end
    
# Elements
Elements_node = aen(Mesh_node,"Elements"; name="Part1", type="hex8")
    for q ∈ eachindex(Em2)
        # aen(Elements_node,"elem",@sprintf("%i, %i, %i, %i, %i, %i, %i, %i",E[q][1],E[q][2],E[q][3],E[q][4],E[q][5],E[q][6],E[q][7],E[q][8]); id = q)
        aen(Elements_node,"elem",join(map(string, Em2[q]), ','); id = q)
    end
Elements_node = aen(Mesh_node,"Elements"; name="Part2", type="hex8")
    for q ∈ eachindex(Em1)
        # aen(Elements_node,"elem",@sprintf("%i, %i, %i, %i, %i, %i, %i, %i",E[q][1],E[q][2],E[q][3],E[q][4],E[q][5],E[q][6],E[q][7],E[q][8]); id = q)
        aen(Elements_node,"elem",join(map(string, Em1[q]), ','); id = q)
    end

    
# Node sets
bcsupport = "bcSuppportList"
bcprescribe = "bcPrescribeList"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ elements2indices(Fb_bottom)],','); name=bcsupport)
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ elements2indices(Fb_top)],','); name=bcprescribe)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    aen(MeshDomains_node,"SolidDomain"; mat = "Normal_material", name="Part1")
    aen(MeshDomains_node,"SolidDomain"; mat = "Multigen_material", name="Part2")
Boundary_node = aen(febio_spec_node, "Boundary")

bc_node = aen(Boundary_node,"bc"; name="FixedDisplacement01", node_set=bcsupport, type="zero displacement")
    aen(bc_node,"x_dof",1)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",1)

LoadData_node = aen(febio_spec_node,"LoadData")

load_controller_node = aen(LoadData_node,"load_controller"; id=1, name="LC_1", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    aen(load_controller_node,"extend","CONSTANT")

points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",0,0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",1,1))

load_controller_node = aen(LoadData_node,"load_controller"; id=2, name="LC_2", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    aen(load_controller_node,"extend","CONSTANT")

points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",1,0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",2,1))

step_node = aen(febio_spec_node,"Step")
step_id = aen(step_node, "step"; id = 1)
control_node = aen(step_id,"Control") # Define Control node: <Control>                  
    aen(control_node,"time_steps",10)
    aen(control_node,"step_size",0.1)

solver_node = aen(control_node,"solver"; type = "solid")
    aen(solver_node,"symmetric_stiffness",1)
    aen(solver_node,"equation_scheme","staggered")
    aen(solver_node,"equation_order","default")
    aen(solver_node,"optimize_bw",0)
    aen(solver_node,"lstol",9e-1)
    aen(solver_node,"lsmin",1e-2)
    aen(solver_node,"lsiter",5)
    # aen(solver_node,"max_refs",70)
    aen(solver_node,"check_zero_diagonal",0)
    aen(solver_node,"zero_diagonal_tol",0)
    aen(solver_node,"force_partition",0)
    aen(solver_node,"reform_each_time_step",1)
    aen(solver_node,"reform_augment",0)
    aen(solver_node,"diverge_reform",1)
    aen(solver_node,"min_residual",1e-20)
    aen(solver_node,"max_residual",0)
    aen(solver_node,"dtol",1e-3)
    aen(solver_node,"etol",1e-2)
    aen(solver_node,"rtol",0)
    aen(solver_node,"rhoi",0)
    aen(solver_node,"alpha",1)
    aen(solver_node,"beta",2.5e-01)
    aen(solver_node,"gamma",5e-01)
    aen(solver_node,"logSolve",0)
    aen(solver_node,"arc_length",0)
    aen(solver_node,"arc_length_scale",0)
qn_method_node = aen(solver_node,"qn_method"; type = "BFGS")
    aen(qn_method_node,"max_ups",0)
    aen(qn_method_node,"max_buffer_size",0)
    aen(qn_method_node,"cycle_buffer",1)
    aen(qn_method_node,"cmax",100000)
time_stepper_node = aen(control_node,"time_stepper"; type = "default")
    aen(time_stepper_node,"max_retries",5)
    aen(time_stepper_node,"opt_iter",6)
    aen(time_stepper_node,"dtmin",1e-3)
    aen(time_stepper_node,"dtmax",0.1)
    aen(time_stepper_node,"aggressiveness",0)
    aen(time_stepper_node,"cutback",5e-1)
    aen(time_stepper_node,"dtforce",0)

aen(control_node,"analysis","STATIC")
aen(control_node,"plot_zero_state",1)
aen(control_node,"plot_range",@sprintf("%.2f, %.2f",0,-1))
aen(control_node,"plot_level","PLOT_MAJOR_ITRS")
aen(control_node,"plot_stride",1)
aen(control_node,"output_level","OUTPUT_MAJOR_ITRS")
aen(control_node,"adaptor_re_solve",1)

boundary = aen(step_id,"Boundary")

bc_name = aen(boundary,"bc"; name = "bcPrescribeList01", node_set = bcprescribe, type = "prescribed displacement")
    aen(bc_name,"dof","z")
    aen(bc_name,"value",displacementMagnitude;lc=@sprintf("%i",1))
    aen(bc_name,"relative",1)

bc_name = aen(boundary,"bc"; name = "FixedDisplacement02", node_set = bcprescribe, type = "zero displacement")
    aen(bc_name,"x_dof",1)
    aen(bc_name,"y_dof",1)
    aen(bc_name,"z_dof",0)

step_id = aen(step_node, "step"; id = 2)
    control_node = aen(step_id,"Control") # Define Control node: <Control>                  
        aen(control_node,"time_steps",10)
        aen(control_node,"step_size",0.1)
    
    solver_node = aen(control_node,"solver"; type = "solid")
        aen(solver_node,"symmetric_stiffness",1)
        aen(solver_node,"equation_scheme","staggered")
        aen(solver_node,"equation_order","default")
        aen(solver_node,"optimize_bw",0)
        aen(solver_node,"lstol",9e-1)
        aen(solver_node,"lsmin",1e-2)
        aen(solver_node,"lsiter",5)
        # aen(solver_node,"max_refs",70)
        aen(solver_node,"check_zero_diagonal",0)
        aen(solver_node,"zero_diagonal_tol",0)
        aen(solver_node,"force_partition",0)
        aen(solver_node,"reform_each_time_step",1)
        aen(solver_node,"reform_augment",0)
        aen(solver_node,"diverge_reform",1)
        aen(solver_node,"min_residual",1e-20)
        aen(solver_node,"max_residual",0)
        aen(solver_node,"dtol",1e-3)
        aen(solver_node,"etol",1e-2)
        aen(solver_node,"rtol",0)
        aen(solver_node,"rhoi",0)
        aen(solver_node,"alpha",1)
        aen(solver_node,"beta",2.5e-01)
        aen(solver_node,"gamma",5e-01)
        aen(solver_node,"logSolve",0)
        aen(solver_node,"arc_length",0)
        aen(solver_node,"arc_length_scale",0)
    qn_method_node = aen(solver_node,"qn_method"; type = "BFGS")
        aen(qn_method_node,"max_ups",0)
        aen(qn_method_node,"max_buffer_size",0)
        aen(qn_method_node,"cycle_buffer",1)
        aen(qn_method_node,"cmax",100000)
    time_stepper_node = aen(control_node,"time_stepper"; type = "default")
        aen(time_stepper_node,"max_retries",5)
        aen(time_stepper_node,"opt_iter",6)
        aen(time_stepper_node,"dtmin",1e-3)
        aen(time_stepper_node,"dtmax",0.1)
        aen(time_stepper_node,"aggressiveness",0)
        aen(time_stepper_node,"cutback",5e-1)
        aen(time_stepper_node,"dtforce",0)
    
    aen(control_node,"analysis","STATIC")
    aen(control_node,"plot_zero_state",1)
    aen(control_node,"plot_range",@sprintf("%.2f, %.2f",0,-1))
    aen(control_node,"plot_level","PLOT_MAJOR_ITRS")
    aen(control_node,"plot_stride",1)
    aen(control_node,"output_level","OUTPUT_MAJOR_ITRS")
    aen(control_node,"adaptor_re_solve",1)
    
    boundary = aen(step_id,"Boundary")
    
    bc_name = aen(boundary,"bc"; name = "bcPrescribeList02", node_set = bcprescribe, type = "prescribed displacement")
        aen(bc_name,"dof","z")
        aen(bc_name,"value",-displacementMagnitude;lc=@sprintf("%i",2))
        aen(bc_name,"relative",1)
    
    bc_name = aen(boundary,"bc"; name = "FixedDisplacement03", node_set = bcprescribe, type = "zero displacement")
        aen(bc_name,"x_dof",1)
        aen(bc_name,"y_dof",1)
        aen(bc_name,"z_dof",0)
            
Output_node = aen(febio_spec_node,"Output")

plotfile_node = aen(Output_node,"plotfile"; type="febio")
    aen(plotfile_node,"var"; type="displacement")
    aen(plotfile_node,"var"; type="stress")
    aen(plotfile_node,"var"; type="relative volume")
    aen(plotfile_node,"var"; type="reaction forces")
    aen(plotfile_node,"var"; type="contact pressure")
    aen(plotfile_node,"compression",@sprintf("%i",0))

logfile_node = aen(Output_node,"logfile"; file=filename_log)
    aen(logfile_node,"node_data"; data="ux;uy;uz", delim=",", file=filename_disp)
    aen(logfile_node,"element_data"; data = "sed", delim =",", file=filename_energy_out)
    # aen(logfile_node,"element_data"; data="s1;s2;s3", delim=",", file=filename_stress)
# <logfile file="tempModel.txt">
#   <node_data data="ux;uy;uz" delim="," file="tempModel_disp_out.txt">1, 2, 3, 4, 5, 6, 7, 8, 

#######

# Write FEB file
XML.write(filename_FEB, doc)

#######
# Run FEBio
run_febio(filename_FEB,FEBIO_PATH)

#######
# Import results
DD = read_logfile(joinpath(saveDir,filename_disp))
DD_stress = read_logfile(joinpath(saveDir,filename_stress))
# DD_E = read_logfile(joinpath(saveDir,filename_energy_out))
#######
##############
filename_disp_new = "febioInputFile_02_DISP.txt"
open(joinpath(saveDir,filename_disp)) do file 
    new_file =[]
    Step_count = 0
    for line in eachline(file)
        if occursin("*Step",line)
            Step_count +=1
            new_line = ("*Step = $Step_count")
            line = new_line 
        end
        push!(new_file,line)
    end 
   
    open(joinpath(saveDir,filename_disp_new),"w") do file_2
        for line in new_file 
            println(file_2,line)
        end
    end

end

filename_E_new = "febioInputFile_energy_out_2.txt"

open(joinpath(saveDir,filename_energy_out)) do file 
    new_file =[]
    Step_count = 0
    for line in eachline(file)
        if occursin("*Step",line)
            Step_count +=1
            new_line = ("*Step = $Step_count")
            line = new_line 
        end
        push!(new_file,line)
    end 
   
    open(joinpath(saveDir,filename_E_new),"w") do file_2
        for line in new_file 
            println(file_2,line)
        end
    end

end
DD_E = read_logfile(joinpath(saveDir,filename_E_new))
# Energy = []
# for i in 1:length(DD_E)
#     E = DD_E[i].data
#     push!(Energy,E)
# end
DD_2 = read_logfile(joinpath(saveDir,filename_disp_new))
time_step = []
for i in 1:length(DD_2)
    time= DD_2[i].time
    push!(time_step,time)
    # println(time)
end
c1_plot = c1.*ones(length(time_step))
cg_plot=c1_g[1].*ones(length(time_step))
cg_plot[time_step .>=1 ] .= c1_g[2]
##################
# Visualization for the boundary
linewidth = 6
strokewidth = 2

fig = Figure(size=(1500,1500))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary Conditions")
hp1=poly!(ax1,GeometryBasics.Mesh(V,Fb), strokewidth=strokewidth,color=:white, strokecolor=:black, shading = FastShading, transparency=true)
hp2 = scatter!(ax1,V[indNodes_bottom],color=:black,markersize=15)
hp3 = scatter!(ax1,V[indNodes_top],color=:red,markersize=15)
Legend(fig[1, 2],[hp2,hp3],["Fixed","Prescribed_disp"])
display(fig)


fig1 = Figure(size=(1500,1500))
  
stepRange = 0:1:length(DD_2)-1
hSlider = Slider(fig1[2, :], range = stepRange, startvalue = length(DD_2)-1,linewidth=30)

nodalColor = lift(hSlider.value) do stepIndex
    norm.(DD_2[stepIndex].data)
end

M = lift(hSlider.value) do stepIndex    
    return GeometryBasics.Mesh(V.+DD_2[stepIndex].data,Fb)
end

titleString = lift(hSlider.value) do stepIndex
  "Step: "*string(stepIndex)
end
C1 = lift(hSlider.value) do stepIndex    
    return c1_plot[stepIndex]
end
Cg = lift(hSlider.value) do stepIndex    
    return cg_plot[stepIndex]
end

Time = lift(hSlider.value) do stepIndex
    return time_step[stepIndex]
end



ax2=Axis3(fig1[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = titleString)
ax3=Axis3(fig1[1, 1], aspect = :data, xlabel = "Time", ylabel = "C_1", title = "Ogden parameter",azimuth=-pi/2,elevation=pi/2)

hp=poly!(ax2,M, color = nodalColor, strokewidth=2, transparency=false, overdraw=false,colormap = Reverse(:Spectral), shading = FastShading)
hp3 = lines!(ax3, time_step, c1_plot, color =:black, linewidth=3)
hp4 = lines!(ax3, time_step, cg_plot, color =:red, linewidth=3)
hp5 = scatter!(ax3, Time,C1, color =:black, markersize = 15)
hp6 = scatter!(ax3, Time,Cg, color =:red, markersize = 15)
Legend(fig1[1, 2],[hp5,hp6],["Materia1","Material2"])

 Colorbar(fig1[1, 4],hp.plots[1],label = "Displacement magnitude [mm]") 

slidercontrol(hSlider,ax2)

display(fig1)

