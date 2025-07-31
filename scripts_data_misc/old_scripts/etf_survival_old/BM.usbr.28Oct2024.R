
#==== 28 Oct 2024 
# Brian Mahardja  ideas

base <- "https://cbr.washington.edu/sac-bin/grow/emergecontrols.pl"

# Describe the output you want. Here: the details of the cohort in a csv file
string <- paste0(base,
                 "?raw=survival",
                 "&reddyear=2012",
                 "&rtnCSV=1")

salmoddefault <- paste0(  "&alevinSurvPar1=2.521",
                          "&alevinSurvPar2=1.461",
                          "&alevinSurvParexp=12",
                          "&eggSurvPar1=1.475",
                          "&eggSurvPar2=1.392",
                          "&eggSurvParexp=11")

waterforumdefault <- paste0("&alevinSurvPar1=1.017554",
                            "&alevinSurvPar2=1.24092",
                            "&alevinSurvParexp=10",
                            "&eggSurvPar1=3.408488",
                            "&eggSurvPar2=1.21122",
                            "&eggSurvParexp=11")


columns <- c(2:11)
morts <- c(paste0("&mortality=exp",salmoddefault),
           paste0("&mortality=exp",waterforumdefault),
           "&mortality=lin")

string2 <- paste0(string,"&units=centigrade&densdeptype=none&tempsource=userdirBrian",
                  "&temfile=ARG2024.csv")
    for(colu in columns){
          for(i in 1:length(morts)){
            mort <- morts[i]
            string3 <- paste0(string2,"&tempdatacolumn=",colu,mort)
            results3 <- scan(string3)
            cat(whichsheet,colu,results3," mort",i,"\n")
          }
    }


