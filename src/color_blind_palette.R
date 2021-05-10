#
# Color-blind friendly colors
#

suppressMessages(require(graphics))


# Palette 1 (8 colors)
# http://jfly.iam.u-tokyo.ac.jp/color/

cbPalette1 <- c("#000000", # 1 Black 
                "#E69F00", # 2 Orange 
                "#56B4E9", # 3 Sky Blue 
                "#009E73", # 4 Bluish Green
                "#F0E442", # 5 Yellow
                "#0072B2", # 6 Blue
                "#D55E00", # 7 Vermillion
                "#CC79A7") # 8 Reddish Purple

palette(cbPalette1)
# plot(1:8,1:8,col=1:8, pch=20,cex=5)
Black_transparency         <- rgb(  0,  0,  0,100,maxColorValue=255)
Orange_transparency        <- rgb(230,159,  0,100,maxColorValue=255)
SkyBlue_transparency       <- rgb( 86,180,233,100,maxColorValue=255)
BluishGreen_transparency   <- rgb(  0,158,115,100,maxColorValue=255)
Yellow_transparency        <- rgb(240,228, 66,100,maxColorValue=255)
Blue_transparency          <- rgb(  0,114,178,100,maxColorValue=255)
Vermillion_transparency    <- rgb(213, 94,  0,100,maxColorValue=255)
ReddishPurple_transparency <- rgb(204,121,167,100,maxColorValue=255)


# Palette 2 (15 colors)

cbPalette2 <- c( rgb(  0,  0,  0, maxColorValue = 255), # 1
                 rgb(  0, 73, 73, maxColorValue = 255), # 2
                 rgb(  0,146,146, maxColorValue = 255), # 3 (Individuals with tritanopia cannot distinguish 3 & 7)
                 rgb(255,109,182, maxColorValue = 255), # 4 (Individuals with tritanopia cannot distinguish 4 & 13)
                 rgb(255,182,119, maxColorValue = 255), # 5
                 rgb( 73,  0,146, maxColorValue = 255), # 6
                 rgb(  0,109,219, maxColorValue = 255), # 7 (Individuals with tritanopia cannot distinguish 3 & 7)
                 rgb(182,109,255, maxColorValue = 255), # 8
                 rgb(109,182,255, maxColorValue = 255), # 9
                 rgb(182,219,255, maxColorValue = 255), # 10
                 rgb(146,  0,  0, maxColorValue = 255), # 11
                 rgb(146, 73,  0, maxColorValue = 255), # 12
                 rgb(219,209,  0, maxColorValue = 255), # 13 (Individuals with tritanopia cannot distinguish 4 & 13)
                 rgb( 36,255, 36, maxColorValue = 255), # 14
                 rgb(255,255,109, maxColorValue = 255)) # 15

# palette(cbPalette2)
# plot(1:15,1:15,col=1:15, pch=20,cex=5)



