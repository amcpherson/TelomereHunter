# Copyright 2015 Lina Sieverling

# This file is part of TelomereHunter.

# TelomereHunter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# TelomereHunter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with TelomereHunter.  If not, see <http://www.gnu.org/licenses/>.



#check if R libraries are installed, otherwise install them
pkgTest <- function(x){
    if (!require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)){
      install.packages(x,dep=TRUE)
      if(!require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)){
        stop(paste("R Package\"", x, "\"not found"))
      } 
    }
  }


for (library in c('ggplot2', 'reshape2', 'grid', 'gridExtra', 'RColorBrewer', 'cowplot', 'svglite')){
  pkgTest(library)
}
