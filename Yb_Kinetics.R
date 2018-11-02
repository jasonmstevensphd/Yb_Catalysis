
library(tidyverse)

# Initial HTE Evaluation of Lewis Acids -----------------------------------

HTE_Data <- read_csv("HTE.csv")

HTE_Data <- HTE_Data %>% rename(`Lewis Acid` = LewisAcid) %>% filter(Entry != "67")

HTE_Data$Product[HTE_Data$Product == "Product 1"] <- 6
HTE_Data$Product[HTE_Data$Product == "Product 2"] <- 7
HTE_Data$Product[HTE_Data$Product == "Product 3"] <- 8

HTE_Plot <- ggplot(HTE_Data, aes(reorder(`Lewis Acid`, `AP PDT`,
                                         FUN = median),
                                 `AP PDT`))+
  geom_boxplot()+
  coord_flip()+
  geom_point(aes(color = Product, shape = `Solvent`), size = 4)+
  labs(x = "Lewis Acid", y = "Area Percent Ester Product")+
  theme(plot.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="grey50"))+
  theme(plot.title = element_text(face = "bold", 
                                  size = 18,
                                  color="Navy"))+
  theme(axis.title = element_text(face = "bold", size = 16))+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust = 1,
                                   size = 12,
                                   angle = 60))+
  theme(axis.text.y = element_text(size = 12))

HTE_Plot

# Ionic Radius Models and Plots -------------------------------------------

Data <- read.csv("Yb_Kinetics.csv")

La_Plot <- Data %>% filter(Lewis_Acid == "La(OTf)3") %>%
  filter(Time_h > 2.5)
Ce_Plot <- Data %>% filter(Lewis_Acid == "Ce(OTf)3") %>%
  filter(Time_h > 2.5)
Pr_Plot <- Data %>% filter(Lewis_Acid == "Pr(OTf)3") %>%
  filter(Time_h > 1.5 & Time_h < 20)
Nd_Plot <- Data %>% filter(Lewis_Acid == "Nd(OTf)3") %>%
  filter(Time_h > 1.5 & Time_h < 20)
Sm_Plot <- Data %>% filter(Lewis_Acid == "Sm(OTf)3") %>%
filter(Time_h > 1.5 & Time_h < 20)
Eu_Plot <- Data %>% filter(Lewis_Acid == "Eu(OTf)3") %>%
  filter(Time_h > 0.5 & Time_h < 8)
Gd_Plot <- Data %>% filter(Lewis_Acid == "Gd(OTf)3") %>%
  filter(Time_h > 0.5 & Time_h < 6)
Tb_Plot <- Data %>% filter(Lewis_Acid == "Tb(OTf)3") %>%
  filter(Time_h > 0.5 & Time_h < 6)
Dy_Plot <- Data %>% filter(Lewis_Acid == "Dy(OTf)3") %>%
  filter(PDT_Yield > 50 & Time_h < 4)
Ho_Plot <- Data %>% filter(Lewis_Acid == "Ho(OTf)3") %>%
  filter(PDT_Yield > 50 & Time_h < 2.5)
Er_Plot <- Data %>% filter(Lewis_Acid == "Er(OTf)3") %>%
  filter(Time_h < 0.75)
Tm_Plot <- Data %>% filter(Lewis_Acid == "Tm(OTf)3") %>%
  filter(Time_h < 0.75)
Yb_Plot <- Data %>% filter(Lewis_Acid == "Yb(OTf)3") %>%
  filter(Time_h < 0.75)
Lu_Plot <- Data %>% filter(Lewis_Acid == "Lu(OTf)3") %>%
  filter(Time_h < 0.75)

Plot <- function(data){
  
  ggplot(data, aes(Time_h, PDT_Yield))+
    geom_point()+
    geom_smooth(method = lm)+
    ggtitle("Product Yield vs Time (h)")+
    labs(x = "Time (h)", y = "In Process Product Yield")+
    theme(plot.title = element_text(face = "bold", 
                                    size = 18,
                                    color="Navy"))+
    theme(axis.title = element_text(face = "bold", size = 16))+
    theme(axis.text.x = element_text(vjust = 1,
                                     hjust = 1,
                                     size = 12,
                                     angle = 60))+
    theme(axis.text.y = element_text(size = 12))
  
}

Model <- function(Data, Yield) {
  
  Model <- lm(Time_h ~ PDT_Yield, Data)
  Model_coeffs <- coefficients(Model)
  Time_to_Completion <- Model_coeffs[1] + Model_coeffs[2]*Yield
  Time_to_Completion

}


Plot(La_Plot)
La_LM <- Model(La_Plot, 95)
La_LM
Data$Time_to_Completion[Data$Lewis_Acid == "La(OTf)3"] <- La_LM

Plot(Ce_Plot)
Ce_LM <- Model(Ce_Plot, 95)
Ce_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Ce(OTf)3"] <- Ce_LM

Plot(Pr_Plot)
Pr_LM <- Model(Pr_Plot, 95)
Pr_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Pr(OTf)3"] <- Pr_LM

Plot(Nd_Plot)
Nd_LM <- Model(Nd_Plot, 95)
Nd_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Nd(OTf)3"] <- Nd_LM

Plot(Sm_Plot)
Sm_LM <- Model(Sm_Plot, 95)
Sm_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Sm(OTf)3"] <- Sm_LM

Plot(Eu_Plot)
Eu_LM <- Model(Eu_Plot, 95)
Eu_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Eu(OTf)3"] <- Eu_LM

Plot(Gd_Plot)
Gd_LM <- Model(Gd_Plot, 95)
Gd_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Gd(OTf)3"] <- Gd_LM

Plot(Tb_Plot)
Tb_LM <- Model(Tb_Plot, 95)
Tb_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Tb(OTf)3"] <- Tb_LM

Plot(Dy_Plot)
Dy_LM <- Model(Dy_Plot, 95)
Dy_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Dy(OTf)3"] <- Dy_LM

Plot(Ho_Plot)
Ho_LM <- Model(Ho_Plot, 95)
Ho_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Ho(OTf)3"] <- Ho_LM

Plot(Er_Plot)
Er_LM <- Model(Er_Plot, 95)
Er_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Er(OTf)3"] <- Er_LM

Plot(Tm_Plot)
Tm_LM <- Model(Tm_Plot, 95)
Tm_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Tm(OTf)3"] <- Tm_LM

Plot(Yb_Plot)
Yb_LM <- Model(Yb_Plot, 95)
Yb_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Yb(OTf)3"] <- Yb_LM

Plot(Lu_Plot)
Lu_LM <- Model(Lu_Plot, 95)
Lu_LM
Data$Time_to_Completion[Data$Lewis_Acid == "Lu(OTf)3"] <- Lu_LM

Data_1 <- Data %>% filter(Lewis_Acid != "Ce(OTf)3" &
                            Lewis_Acid != "Y(OTf)3") %>%
  mutate(log_Time_to_Completion = log10(Time_to_Completion),
         log_Ionic_Radius = log10(Ionic_Radius))

Ionic_Radius_Model <- lm(Ionic_Radius ~ Time_to_Completion, Data_1)
summary(Ionic_Radius_Model)


Data_1$Annotate[Data_1$Lewis_Acid == "La(OTf)3"] <- "La"
Data_1$Annotate[Data_1$Lewis_Acid == "Pr(OTf)3"] <- "Pr"
Data_1$Annotate[Data_1$Lewis_Acid == "Nd(OTf)3"] <- "Nd"
Data_1$Annotate[Data_1$Lewis_Acid == "Sm(OTf)3"] <- "Sm"
Data_1$Annotate[Data_1$Lewis_Acid == "Eu(OTf)3"] <- "Eu"
Data_1$Annotate[Data_1$Lewis_Acid == "Gd(OTf)3"] <- "Gd"
Data_1$Annotate[Data_1$Lewis_Acid == "Tb(OTf)3"] <- "Tb"
Data_1$Annotate[Data_1$Lewis_Acid == "Dy(OTf)3"] <- "Dy"
Data_1$Annotate[Data_1$Lewis_Acid == "Ho(OTf)3"] <- "Ho"
Data_1$Annotate[Data_1$Lewis_Acid == "Er(OTf)3"] <- "Er"
Data_1$Annotate[Data_1$Lewis_Acid == "Tm(OTf)3"] <- "Tm"
Data_1$Annotate[Data_1$Lewis_Acid == "Yb(OTf)3"] <- "Yb"
Data_1$Annotate[Data_1$Lewis_Acid == "Lu(OTf)3"] <- "Lu"


Ionic_Radius_Plot <- ggplot(Data_1, aes(Ionic_Radius, Time_to_Completion))+
  geom_point(size = 2)+
  geom_text(aes(label = Annotate),
             nudge_y = 2)+
  ggtitle("Ionic Radius vs Time to Completion (h)")+
  labs(x = "Ionic Radius of Lanthanide(III) Triflate (Å)", y = "Time to Reaction Completion (h)")+
  theme(plot.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="grey50"))+
  theme(plot.title = element_text(face = "bold", 
                                  size = 18,
                                  color="Navy"))+
  theme(axis.title = element_text(face = "bold", size = 16))+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust = 1,
                                   size = 12,
                                   angle = 60))+
  theme(axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = round(seq(min(Data_1$Ionic_Radius), max(Data_1$Ionic_Radius), by = 1),1))+
  scale_y_continuous(breaks = round(seq(min(Data_1$Time_to_Completion), max(Data_1$Time_to_Completion), by = 2.5),2.5))

Ionic_Radius_Plot

Ionic_Radius_Plot_1 <- ggplot(Data_1, aes(log_Ionic_Radius, log_Time_to_Completion))+
  geom_point()+
  ggtitle("log(Ionic Radius) vs log(Time to Completion)")+
  labs(x = "log(Ionic Radius) (Å)", y = "log(Time to Completion (h))")+
  theme(plot.title = element_text(face = "bold", 
                                  size = 18,
                                  color="Navy"))+
  theme(axis.title = element_text(face = "bold", size = 16))+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust = 1,
                                   size = 12,
                                   angle = 60))+
  theme(axis.text.y = element_text(size = 12))

Ionic_Radius_Plot_1


# Kinetic Isotope Effect --------------------------------------------------

library(readxl)

KIE_Data <- read_xlsx("Yb_KIE.xlsx", sheet = "Sheet1")

KIE_Data <- KIE_Data %>% filter(Time_h <3)

KIE_Data <- KIE_Data %>% rename(`Product (S)-6` = Product_3) %>%
  rename(`Product 23` = Product_25)

KIE_Plot <- ggplot(KIE_Data, aes())+
  geom_point(aes(x = Time_h, y = `Product (S)-6`, color = "Product (S)-6"), size = 2)+
  geom_smooth(aes(x = Time_h, y = `Product (S)-6`), method = "lm", color = "red")+
  geom_point(aes(x = Time_h, y = `Product 23`, color = "Product 23"), size = 2)+
  geom_smooth(aes(x = Time_h, y = `Product 23`), method = "lm", color = "blue")+
  scale_colour_manual(name="legend", values=c("red", "blue"))+
  ggtitle("Reaction Progress for Methanol\nvs Methanol-d4")+
  labs(x = "Time (h)", y = "Area Percent Product")+
  theme(plot.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="grey50"))+
  theme(plot.title = element_text(face = "bold", 
                                  size = 18,
                                  color="Navy"))+
  theme(axis.title = element_text(face = "bold", size = 16))+
  theme(axis.text.x = element_text(vjust = 1,
                                   hjust = 1,
                                   size=12))+
  theme(axis.text.y = element_text(size=12))

KIE_Plot

sessionInfo()
