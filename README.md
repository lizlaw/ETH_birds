# ETH_birds
Multi-species occupancy model of birds in Jimma zone, ETH
Contact: Liz Law workingconservation@gmail.com

_Objective_: to model species richness of birds across Jimma zone, ETH. Data include 1ha point observations for 80 farm and 65 forest sites, each with two visits. Total of 159 species (observations observed to genus only considered as additional species, identifiable as _ spp in codes). Models will be projected to the entire study region, under current and future land use scenarios. 

_Methods_: Models developed as multi-species occupancy models (in nimble), as these are able to account for imperfect observation, and leverage data from common species to inform rarer ones. Individual-based rarefication for each habitat type (code elsewhere) was used to set the maximum number of species observable in each habitat type. Prior cleaning removed several sites for which only one visit was available, or for which the habitat (farm/forest) did not match spatial variables. 

Initial variables considered include:
  * Occupancy: elevation, topographic wetness index, heat load index, slope, distance to forest edge, historical distance to forest edge, farm/forest type (i.e. whether farm/forest in 1985), percent woody vegetation, landscape diversity (simpsons diversity) at 1ha and 200m diameter.  
  * Detection: start time (including quadratic), date, visibility, cloud, number of observers, recording.
These initial variables were transformed/re-categorised and scaled+centred to improve normality/spread, then pre-selected via pearson's correlation, stepwise AIC glm, random forest, and multi-model inference (on glm models). Starting models for multi-occupancy:
  * FARM: elevation+fl_dis+sidi1ha+slope+farm_type || date+start1+start2+visibility
  * FOREST: elevation+fr_dis+slope+forest_type || date+start1+start2+n_observers

Poor convergence of initial models likely due to overparameterisation, and currently attempting to simplify the model. Versions:

  * v5. initial model versions with extraneous predictor variables removed
  * v6. removed tau (not required in nimble), omega (doesn't do anything here, as we don't select from a regional pool - we estimate directly the number of missing species via rarefication), assume that all npsi (occ parameters) share the same sd.betalpsi, assume that all np (obs parameters) share the same sd.betalp
  * v7. v6 + assumed that all rare species (observed <1 in all study) have the same same beta parameters. 

All code contained within ~R/ folder

Note: gitignore set for ~Data/ and ~Results/ folders because of data permissions and large data samples produced
