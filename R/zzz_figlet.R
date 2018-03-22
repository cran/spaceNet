# .onAttach <- function(lib, pkg)
# {
#   version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
#   if(interactive())
#   { # colossal -- try also Nancyj
#     packageStartupMessage(
# "                                            888b    888          888
#                                             8888b   888          888
#                                             88888b  888          888
# .d8888b  88888b.   8888b.   .d8888b .d88b.  888Y88b 888  .d88b.  888888
# 88K      888 \"88b     \"88b d88P\"   d8P  Y8b 888 Y88b888 d8P  Y8b 888
# \"Y8888b. 888  888 .d888888 888     88888888 888  Y88888 88888888 888
#      X88 888 d88P 888  888 Y88b.   Y8b.     888   Y8888 Y8b.     Y88b.
#  88888P\' 88888P\"  \"Y888888  \"Y8888P \"Y8888  888    Y888  \"Y8888   \"Y888
#          888
#          888
#          888                                                            ",
#       "version ", version, "\n" )
#   }
#   else
#   { packageStartupMessage("Package 'spaceNet' version ", version) }
#
#   packageStartupMessage("Type 'citation(\"spaceNet\")' for citing this R package in publications.")
#   invisible()
# }


.onAttach <- function(lib, pkg)
{
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  if(interactive())
  { # colossal -- try also Nancyj
    packageStartupMessage(
      "
                                             888    8b            88
.d8888b. 88d888b. .d8888b. .d8888b. .d8888b. 88 8   88 .d8888b. d8888P
Y8ooooo. 88'  `88 88'  `88 88       88ooood8 88  8  88 88ooood8   88
      88 88.  .88 88.  .88 88.  ... 88.  ... 88   8 88 88.  ...   88
`88888P' 88Y888P' `88888P8 `88888P   88888P  dP     dP  88888P    dP
         88
         dP                                                            ",
      "version ", version, "\n" )
  }
  else
  { packageStartupMessage("Package 'spaceNet' version ", version) }

  packageStartupMessage("Type 'citation(\"spaceNet\")' for citing this R package in publications.")
  invisible()
}
