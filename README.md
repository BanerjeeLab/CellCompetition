# CellCompetition
Cell-based model for cell competition.

These scripts run on the program CompuCell3D version 4.5 (downloadable from https://compucell3d.org/SrcBin).
CompuCell3D is an open source Cellular Potts Model simulator maintained by Prof. James Glazier (& friends). 
Instructions and tutorials on how to run CompuCell3D are provided on their website https://compucell3d.org


### **An overvie of the basic content in each file** <br/>
Alter the rules and parameters of the simulations in:
<pre> Competition_Scripts\Simulation\HST_Steppables.py </pre>

Choose the parameter search space in:
<pre> Competition_Scripts\Simulation\ParameterScanSpecs.json </pre>

Choose which phenomena you model by changing which classes you call in:
<pre> Competition_Scripts\Simulation\Comp_Conf.py </pre>

Define simulation parameters in :
<pre> Competition_Scripts\Simulation\Comp_Conf.xml </pre>

---------------------------------------------------------------------------------------

Before running these scripts you will need to unzip the piff file located here:  <br/>  <pre> Competition_Scripts\Simulation\init_conf.zip </pre> <br/>
The final unzipped file must be located here:  <br/> <pre> Competition_Scripts\Simulation\init_conf.piff </pre>
