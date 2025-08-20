from cc3d import CompuCellSetup



from HST_Steppables import Initializer
CompuCellSetup.register_steppable( steppable = Initializer(frequency=1, piff_start=True))

from HST_Steppables import Simulation_Killer
CompuCellSetup.register_steppable( steppable = Simulation_Killer(frequency=1))

from HST_Steppables import Crowded_Growth_Steppable
CompuCellSetup.register_steppable( steppable = Crowded_Growth_Steppable(frequency=1))



from HST_Steppables import Sizer_Timer_Mitosis_Steppable
CompuCellSetup.register_steppable( steppable = Sizer_Timer_Mitosis_Steppable(frequency=10))

from HST_Steppables import Death_Steppable
CompuCellSetup.register_steppable( steppable = Death_Steppable(frequency=10))

from HST_Steppables import Plot_Total_Div_and_Apo_Counts
CompuCellSetup.register_steppable( steppable = Plot_Total_Div_and_Apo_Counts(frequency=10))

from HST_Steppables import Extrusion
CompuCellSetup.register_steppable(steppable=Extrusion(frequency=10))



# from HST_Steppables import RadialForce
# RadialForceInstance = RadialForce(sim,_frequency=5)
# CompuCellSetup.registerSteppable( RadialForceInstance )

# from HST_Steppables import ColorCellsByCellCyclePhaseSteppable
# ColorPhaseInstance = ColorCellsByCellCyclePhaseSteppable(sim,_frequency=10)
# CompuCellSetup.registerSteppable( ColorPhaseInstance )

# from HST_Steppables import ColorCellsByGrowthRateSteppable
# ColorGInstance = ColorCellsByGrowthRateSteppable(sim,_frequency=1)
# CompuCellSetup.registerSteppable( ColorGInstance )

# from HST_Steppables import Color_Cells_Switch_Between_GR_and_Cycle_Phase_Steppable
# ColorGInstance = Color_Cells_Switch_Between_GR_and_Cycle_Phase_Steppable(sim,_frequency=1, _kill_sim_after_colony_area_threshold_met=False)
# CompuCellSetup.registerSteppable( ColorGInstance )

CompuCellSetup.run()        