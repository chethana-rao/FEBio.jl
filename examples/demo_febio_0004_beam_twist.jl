using Comodo # https://github.com/COMODO-research/Comodo.jl
using FEBio
using Printf
using GLMakie
using LinearAlgebra
using XML
using GeometryBasics
using Statistics

######
# Set FEBio path here for now
const FEBIO_PATH = "/home/simulimb/FEBioStudio/bin/febio4" # Path to FEBio executable

###### 
# Control parameters 
pointSpacing = 3.0
rotation_prescribed = 6.2831853e+00

# Material parameters
c1 = 1e-3 #Shear-modulus-like parameter
m1 = 2  #Material parameter setting degree of non-linearity
k_factor = 1e2  #Bulk modulus factor 
k = c1*k_factor #Bulk modulus

###### 
# Creating a hexahedral mesh for a cube 
boxDim = [10.0,40.0,10.0] # Dimensionsions for the box in each direction
boxEl = ceil.(Int64,boxDim./pointSpacing) # Number of elements to use in each direction 
E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl)

# Create face sets to define node sets later 
Fb_bottom = Fb[CFb_type.==1]

Fb_s1 = Fb[CFb_type.==6]
Fb_front = Fb[CFb_type.==3]

indNodesBack = elements2indices(Fb[CFb_type.==4])
indNodesFront = elements2indices(Fb[CFb_type.==3])
# indRigid = elements2indices(Fb_bottom)
COM = mean(V[indNodesFront])
######
# Define file names
saveDir = joinpath(febiojl_dir(),"assets","temp") # Main directory to save FEBio input and output files
if !isdir(saveDir)
    mkdir(saveDir)      
end

filename_FEB = joinpath(saveDir,"febioInputFile_01.feb")   # The main FEBio input file
filename_xplt = joinpath(saveDir,"febioInputFile_01.xplt") # The XPLT file for viewing results in FEBioStudio
filename_log = joinpath(saveDir,"febioInputFile_01_LOG.txt") # The log file featuring the full FEBio terminal output stream
filename_disp = "febioInputFile_01_DISP.txt" # A log file for results saved in same directory as .feb file  e.g. nodal displacements
filename_stress = "febioInputFile_01_STRESS.txt"

######
# Define febio input file XML
doc,febio_spec_node = feb_doc_initialize()

aen(febio_spec_node,"Module"; type = "solid") # Define Module node: <Module type="solid"/>

control_node = aen(febio_spec_node,"Control") # Define Control node: <Control>
    aen(control_node,"analysis","STATIC")               
    aen(control_node,"time_steps",30)
    aen(control_node,"step_size",3.33e-2)
    aen(control_node,"plot_zero_state",1)
    aen(control_node,"plot_range",@sprintf("%.2f, %.2f",0,-1))
    aen(control_node,"plot_level","PLOT_MAJOR_ITRS")
    aen(control_node,"plot_stride",1)
    aen(control_node,"output_level","OUTPUT_MAJOR_ITRS")
    aen(control_node,"adaptor_re_solve",1)

time_stepper_node = aen(control_node,"time_stepper"; type = "default")
    aen(time_stepper_node,"max_retries",5)
    aen(time_stepper_node,"opt_iter",10)
    aen(time_stepper_node,"dtmin",1e-3)
    aen(time_stepper_node,"dtmax",0.1)
    aen(time_stepper_node,"aggressiveness",0)
    aen(time_stepper_node,"cutback",5e-1)
    aen(time_stepper_node,"dtforce",0)

solver_node = aen(control_node,"solver"; type = "solid")
    aen(solver_node,"symmetric_stiffness",1)
    aen(solver_node,"equation_scheme",1)
    aen(solver_node,"equation_order","default")
    aen(solver_node,"optimize_bw",0)
    aen(solver_node,"lstol",9e-1)
    aen(solver_node,"lsmin",1e-2)
    aen(solver_node,"lsiter",5)
    aen(solver_node,"max_refs",25)
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

Globals_node   = aen(febio_spec_node,"Globals")

Constants_node = aen(Globals_node,"Constants")
    aen(Constants_node,"R",8.3140000e-06)
    aen(Constants_node,"T",298)
    aen(Constants_node,"F",9.6485000e-05)

Material_node = aen(febio_spec_node,"Material")

material_node = aen(Material_node,"material"; id = 1, name="Material1", type="Ogden")
    aen(material_node,"c1",c1)
    aen(material_node,"m1",m1)
    aen(material_node,"c2",c1)
    aen(material_node,"m2",-m1)
    aen(material_node,"k",k)

material_node = aen(Material_node,"material"; id = 2, name="Material2", type="rigid body")
    aen(material_node,"density",1)
    aen(material_node,"center_of_mass",join([@sprintf("%e, %e, %e",COM[1],COM[2],COM[3])]))
   
Mesh_node = aen(febio_spec_node,"Mesh")

Nodes_node = aen(Mesh_node,"Nodes"; name="nodeSet_all")
    for q ∈ eachindex(V)
        # aen(Nodes_node,"node",@sprintf("%.2f, %.2f, %.2f",V[q][1],V[q][2],V[q][3]); id = q)
        aen(Nodes_node,"node", join([@sprintf("%.16e",x) for x ∈ V[q]],','); id = q)     
    end
    
# Elements
Elements_node = aen(Mesh_node,"Elements"; name="Part1", type="hex8")
    for q ∈ eachindex(E)
        # aen(Elements_node,"elem",@sprintf("%i, %i, %i, %i, %i, %i, %i, %i",E[q][1],E[q][2],E[q][3],E[q][4],E[q][5],E[q][6],E[q][7],E[q][8]); id = q)
        aen(Elements_node,"elem",join(map(string, E[q]), ','); id = q)
    end
Elements_node = aen(Mesh_node,"Elements"; name="Part2", type="quad4")
for q ∈ eachindex(Fb_front)
    # aen(Elements_node,"elem",@sprintf("%i, %i, %i, %i, %i, %i, %i, %i",E[q][1],E[q][2],E[q][3],E[q][4],E[q][5],E[q][6],E[q][7],E[q][8]); id = q)
    aen(Elements_node,"elem",join(map(string, Fb_front[q]), ','); id = q+size(E)[1])
end
    
# Node sets
bcSupportList = "bcSupportList"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ indNodesBack],','); name=bcSupportList)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")
    aen(MeshDomains_node,"ShellDomain"; mat = "Material2", name="Part2")
Boundary_node = aen(febio_spec_node, "Boundary")

bc_node = aen(Boundary_node,"bc"; name="zero_displacement_xyz", node_set=bcSupportList, type="zero displacement")
    aen(bc_node,"x_dof",1)
    aen(bc_node,"y_dof",1)
    aen(bc_node,"z_dof",1)

Rigid_node = aen(febio_spec_node, "Rigid")

rig_node = aen(Rigid_node,"rigid_bc"; name= "RotFix_1", type ="rigid_fixed") 
    aen(rig_node,"rb",2)
    aen(rig_node,"Rx_dof",1)
    aen(rig_node,"Ry_dof",1)
    aen(rig_node,"Rz_dof",1)
    aen(rig_node,"Ru_dof",1)
    aen(rig_node,"Rw_dof",1)
rig_node = aen(Rigid_node,"rigid_bc"; name = "RigidPrescribe", type ="rigid_rotation")
    aen(rig_node,"rb",2)
    aen(rig_node,"dof","Rv")
    aen(rig_node,"value",rotation_prescribed; lc=@sprintf("%i",1))    
    aen(rig_node,"relative",0)


LoadData_node = aen(febio_spec_node,"LoadData")

load_controller_node = aen(LoadData_node,"load_controller"; id=1, name="LC_1", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    
points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",0,0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",1,1))

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
    aen(logfile_node,"element_data"; data="s1;s2;s3", delim=",", file=filename_stress)
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

#######
# Visualization
linewidth = 6
strokewidth = 2
markersize = 20

fig = Figure(size=(800,800))

ax1=Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary Conditions")
hp1=poly!(ax1,GeometryBasics.Mesh(V,Fb), strokewidth=strokewidth,color=:white, strokecolor=:black, shading = FastShading, transparency=true)
hp2 = scatter!(ax1,V[indNodesBack],color=:black,markersize=markersize)
hp3 = scatter!(ax1,V[indNodesFront],color=:red,markersize=markersize)
Legend(fig[1, 2],[hp2,hp3],["Fixed","Rigid"])

stepRange = 0:1:length(DD)-1
hSlider = Slider(fig[2, 3], range = stepRange, startvalue = length(DD)-1,linewidth=30)

nodalColor = lift(hSlider.value) do stepIndex
    norm.(DD[stepIndex].data)
end

M = lift(hSlider.value) do stepIndex    
    return GeometryBasics.Mesh(V.+DD[stepIndex].data,Fb)
end

titleString = lift(hSlider.value) do stepIndex
  "Step: "*string(stepIndex)
end

ax2=Axis3(fig[1, 3], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Deformed State: $titleString")

hp=poly!(ax2,M, strokewidth=2,color=nodalColor, transparency=false, overdraw=false,colormap = Reverse(:Spectral), shading = FastShading)
Colorbar(fig[1, 4],hp.plots[1],label = "Displacement magnitude [mm]") 

slidercontrol(hSlider,ax2)

fig

