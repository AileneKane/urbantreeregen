## README for urbantreeregen repo
## Ailene Ettinger, ailene.ettinger@tnc.org
#############################################
# Contains data and code associated with the following reference:
# Ettinger, A.K., Lee, B.R. and Montgomery, S., 2017. Seed limitation and lack of downed wood, not invasive species, threaten conifer regeneration in an urban forest. Urban ecosystems, 20, pp.877-887.
# Available at https://link.springer.com/article/10.1007/s11252-016-0640-3
# For all datafiles, we use the following abbreviates:
# 1) For "treatment" column:NI=no ivy (ivy removed), Control = ivy present (not removed) and no wood chips added, WCNI=woodchips added, no ivy (ivy removed), WC= wood chips added, ivy present (not removed)
# 2) For "species" column: ABGR= Abies grandis, PSME=PSeudotsuga menziesii, THPL = Thuja plicata, TSHE= Tsuga heterophylla
# Data files include the following:
# 1) LPGERM20112012.csv contains germination data
###### "numseeds" column = number of seeds planted
###### "numgerm" column = number of seeds that germinated at some point during the study period
###### "totalfails" column = number of seeds that did not germinate (totalfails=numseeds-numgerm)
###### "survivinggerm" column = number of seeds that had transitioned to living seedlings and were still alive at the end of the year 
###### "failedgerm" column = number of seeds that germinated but were not alive at the end of the year (failedgerm = numgerm-survivinggerm)

# 2)  TSHEtransplant.csv and THPLtransplant.csv contains survival and height data for TSHE and THPL seedlings, respectively, from 2011-2012
###### For "Wood" and "Ivy" columns, 0=absent (wood chips not added, or ivy removed) and 1=present  (woodchips added, or ivy not removed)
###### The "SeedlingNo" column is a unique id/tag number and color assigned to the seedling
###### For all "Ht" columns, the date is contained in the title, units are centimeters
###### For all "Status" columns, 0=dead, 1=alive
###### For Notes column, "TR"= tag removed

# 3) LPMORAseedtrap.csv contains seed trap data from Lincoln park, Seattle, WA (LP), the site of the transplant study and from three old growth sites at Mount Rainier National Park (TA01, TO11, TO04)



