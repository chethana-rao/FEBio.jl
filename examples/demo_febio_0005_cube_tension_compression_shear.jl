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
pointSpacing = 2.0
stretchLoad = 1.3
displacementMagnitude = (stretchLoad*sampleSize) - sampleSize 

#Material Constants
c1 = 1e-3 #Shear-modulus-like parameter
m1 = 8  #Material parameter setting degree of non-linearity
k_factor = 1e2  #Bulk modulus factor 
k = c1*k_factor #Bulk modulus

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
    
# Node sets
bcsupport = "bcSuppportList"
bcprescribe = "bcPrescribeList"
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ elements2indices(Fb_bottom)],','); name=bcsupport)
aen(Mesh_node,"NodeSet",join([@sprintf("%i",x) for x ∈ elements2indices(Fb_top)],','); name=bcprescribe)

MeshDomains_node = aen(febio_spec_node, "MeshDomains")
    aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")

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

load_controller_node = aen(LoadData_node,"load_controller"; id=3, name="LC_3", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    aen(load_controller_node,"extend","CONSTANT")

points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",2,0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",3,1))
load_controller_node = aen(LoadData_node,"load_controller"; id=4, name="LC_4", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    aen(load_controller_node,"extend","CONSTANT")
points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",3,0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",4,1))

load_controller_node = aen(LoadData_node,"load_controller"; id=5, name="LC_5", type="loadcurve")
    aen(load_controller_node,"interpolate","LINEAR")
    aen(load_controller_node,"extend","CONSTANT")
points_node = aen(load_controller_node,"points")
    aen(points_node,"pt",@sprintf("%.2f, %.2f",4,0))
    aen(points_node,"pt",@sprintf("%.2f, %.2f",5,1))

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

bc_name = aen(boundary,"bc"; name = "bcPrescribeList", node_set = bcprescribe, type = "prescribed displacement")
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
    
    bc_name = aen(boundary,"bc"; name = "bcPrescribeList", node_set = bcprescribe, type = "prescribed displacement")
        aen(bc_name,"dof","z")
        aen(bc_name,"value",-displacementMagnitude;lc=@sprintf("%i",2))
        aen(bc_name,"relative",1)
    
    bc_name = aen(boundary,"bc"; name = "FixedDisplacement03", node_set = bcprescribe, type = "zero displacement")
        aen(bc_name,"x_dof",1)
        aen(bc_name,"y_dof",1)
        aen(bc_name,"z_dof",0)
step_id = aen(step_node, "step"; id = 3)
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

bc_name = aen(boundary,"bc"; name = "bcPrescribeList", node_set = bcprescribe, type = "prescribed displacement")
    aen(bc_name,"dof","z")
    aen(bc_name,"value",-displacementMagnitude;lc=@sprintf("%i",3))
    aen(bc_name,"relative",1)

bc_name = aen(boundary,"bc"; name = "FixedDisplacement04", node_set = bcprescribe, type = "zero displacement")
    aen(bc_name,"x_dof",1)
    aen(bc_name,"y_dof",1)
    aen(bc_name,"z_dof",0)

step_id = aen(step_node, "step"; id = 4)
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
    
    bc_name = aen(boundary,"bc"; name = "bcPrescribeList", node_set = bcprescribe, type = "prescribed displacement")
        aen(bc_name,"dof","z")
        aen(bc_name,"value",displacementMagnitude;lc=@sprintf("%i",4))
        aen(bc_name,"relative",1)
    
    bc_name = aen(boundary,"bc"; name = "FixedDisplacement05", node_set = bcprescribe, type = "zero displacement")
        aen(bc_name,"x_dof",1)
        aen(bc_name,"y_dof",1)
        aen(bc_name,"z_dof",0)

step_id = aen(step_node, "step"; id = 5)
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
        
        bc_name = aen(boundary,"bc"; name = "bcPrescribeList", node_set = bcprescribe, type = "prescribed displacement")
            aen(bc_name,"dof","z")
            aen(bc_name,"value",displacementMagnitude;lc=@sprintf("%i",5))
            aen(bc_name,"relative",1)
        
        bc_name = aen(boundary,"bc"; name = "FixedDisplacement06", node_set = bcprescribe, type = "zero displacement")
            aen(bc_name,"x_dof",0)
            aen(bc_name,"y_dof",0)
            aen(bc_name,"z_dof",1)
            
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
hp2 = scatter!(ax1,V[indNodes_bottom],color=:black,markersize=markersize)
hp3 = scatter!(ax1,V[indNodes_top],color=:red,markersize=markersize)
Legend(fig[1, 2],[hp2,hp3],["Fixed","Prescribed_disp"])



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

