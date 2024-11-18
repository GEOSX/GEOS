#!/usr/bin/perl
use strict;

my $narg = $#ARGV + 1;

die "usage: $0 <use_y>\n" if($narg != 1);

my $useY = shift;


my $dimMinSet = 5;
my $dimMaxSet = 4;
if($useY==1)
{
    my $dimMinSet += 2;
    my $dimMaxSet += 2;    
}


print q{<?xml version="1.0" ?>
<!-- GPAC.x -i mesh.xml -->
<Problem>

  <Units length="m" mass="kg" time="us" temperature="K" mole="mol" />

  <!-- MESH -->
  <Parameter meshFile="mesh.geom"/>
  <Mesh file="$:meshFile" units="cm" />
  <ElementRegions>
    <EB1 elementtype="uniformstrain" BulkModulus="16.0 GPa" ShearModulus="12.0 GPa" Density="2000" hgDamp="0.1" hgStiff="0.01"/>
    <EB2 elementtype="uniformstrain" BulkModulus="16.0 GPa" ShearModulus="12.0 GPa" Density="2000" hgDamp="0.1" hgStiff="0.01"/>
  </ElementRegions>

  <!-- CONTACT -->
  <Contact active="1" 
	   aperture="1.0e-2 mm" 
	   apertureFactor="0.1"
	   maximumAperture="0.2 mm"
	   arealStiffnessMin="4.0 MPa"
	   arealStiffnessMax="0.4 GPa"
           arealShearStiffness="4.0 MPa"
	   cosMinTol="0.7"
   />

  <!-- SOLVERS -->
  <Solvers>

    <LagrangeExplicitDynamicsSolver name="solver1" courant="0.5" />

    <SteadyStateParallelPlateFlowSolver_TwoD name="pvSolver"
                                             tol="1e-10" 
                                             faceset="NS2" verbose="0" 
                                             MinimumAperture="1e-7 mm"
                                             MaximumAperture="0.2 mm"/>
  </Solvers>

  <!-- SOLVER APPLICATIONS -->
  <SolverApplications>
    <Application name="1" begintime="0.0" endtime="1e4 us">
      <Apply solver="solver1" toregions="EB1 EB2"/>
      <Apply solver="pvSolver" toregions="EB2" />
    </Application>  
  </SolverApplications>

  <!-- BOUNDARY CONDITIONS -->
  <BoundaryConditions>

    <!-- Fracture flow boundary conditions -->
    <BoundaryCondition     object="Edge"
                           fieldname="Pressure"
}; 
print "                           setname=\"NS$dimMinSet\"\n";
print q{                           scale="1.0 MPa" />
    <BoundaryCondition     object="Edge"
                           fieldname="Pressure"
};
print "                           setname=\"NS$dimMaxSet\"\n";
print q{                           scale="0.0 MPa" />

    <!-- Mechanical boundary conditions -->
    <BoundaryCondition  fieldname="Velocity"
                        setname="NS1"
                        component="0" 
                        scale="-1.0"
                        timetable="tv1" />    
    <BoundaryCondition  fieldname="Velocity"
                        setname="NS2"
                        component="0" 
                        scale="1.0"
                        timetable="tv1" />
    <TractionBoundaryCondition setname="NS1" 
	                       direction="0.0 0.0 -1.0" 
			       scale="2.0e6"
			       timetable="tractionTable" />
    <TractionBoundaryCondition setname="NS2" 
	                       direction="0.0 0.0 1.0" 
			       scale="2.0e6"
			       timetable="tractionTable" />
  </BoundaryConditions>

  <!-- INITIAL CONDITIONS-->
  <InitialConditions>
    <ConstantInitialCondition setname="NS3" fieldname="FlowSurface" object="Face" value="1.0" />
    <ConstantInitialCondition fieldname="Pressure" fieldtype="Scalar" object="Face" value="1.0" />

    <!-- Aperture -->
    <CalculateFaceCenters/>
    <CalculateAperture/> 

  </InitialConditions>

  <!-- TABLES -->
  <Tables>
    <Table1D name="tractionTable"
               coord  = "0.0, 0.001,  1.0e9"
               value = "0.0,    1.0,    1.0" />
    <Table1D name="tv1"
               coord  = "0.0, 0.001,  0.00101, 1.0e9"
               value = "0.0,    0.0,    1.0,   1.0" />
    <Table1D name="ts"
               coord  = "0.0, 0.00001,  1.0e9"
               value = "1.0,    1.0,    1.0" />
  </Tables>

  <!-- SPATIAL DECOMPOSITION-->
  <Partition>
    <SpatialPartition xpar="2" ypar="1" zpar="1" />
  </Partition>
    

  <!-- OUTPUT -->
  <Output  plotfile_root="silo_" parallel_silo="1" plot_interval="0.079 us" writeFEMEdges="1" writeFEMFaces="1"/>

</Problem>
};
