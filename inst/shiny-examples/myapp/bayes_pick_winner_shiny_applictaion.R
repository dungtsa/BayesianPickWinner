
library(shiny)
require(clinfun)
library(knitr)
library(rmarkdown)


Bayesian_posterior_probability <-function(n_response_armA, n_nonresponse_armA, n_response_armB, n_nonresponse_armB,sim.n=100000)
{
  # number of response in Arm A: n_response_armA
  # number of non-response in Arm A: n_nonresponse_armA
  # number of response in Arm B: n_response_armB
  # number of non-response in Arm B: n_nonresponse_armB
  # n.sim: the number of simulations

  response_rate_A <-rbeta(sim.n,  1 + n_response_armA,
                          1 + n_nonresponse_armA);

  response_rate_B <- rbeta(sim.n, 1 + n_response_armB,
                           1 + n_nonresponse_armB)

  mean(response_rate_B > response_rate_A)  # ---- this is Bayesian posterior probability, Pr(B>A)

}




Bayesian.simon.two.stage.fun<-function(sim.n0=10, study.design='optimal',  beta1=.1,alpha1=.1, p10=0.2,p20=0.4,n.max=NULL)
{
  require(clinfun)
#  set.seed(1000)

  if(study.design=='optimal') optimal<-T else optimal<-F

  prob.mean.v<-function(data)
  {
    sim.n<-100000
    aa1<-data[1] # number of response in Arm A
    bb1<-data[2] # number of non-response in Arm A
    aa2<-data[3] # number of response in Arm B
    bb2<-data[4] # number of non-response in Arm B
    beta1<-rbeta(sim.n,aa1,bb1);
    beta2<-rbeta(sim.n,aa2,bb2)
    mean(beta1<beta2)
  }

  fun.bay.fisher.v<-function(x,prior_armA=c(1,1),prior_armB=c(1,1)){
    table1<-cbind(x[1:2],x[3:4])
    # x[1:2]= (# of response, # of non-response) in arm A
    # x[3:4]= (# of response, # of non-response) in arm B
    # prior_armA= two beta parameters: a and b for Arm A
    # prior_armB= two beta parameters: a and b for Arm B

    p.fisher.less<-fisher.test(table1,alternative ='less')$p.value
    c(p.fisher.less,prob.mean.v(x+c(prior_armA,prior_armB)))
  }

  baye.sim.fun<-function(sim.n0,p1.true,p2.true,n1, n.total,n.stop1,n.stop2)
  {
    ans<-numeric(0)

    for(i in 1:sim.n0)
    {
      x1<-rbinom(n.total,1,p1.true)
      x2<-rbinom(n.total,1,p2.true)

      ans.A<-ans.B<-NULL
      p.tmp<-poster.mean1<-NA
      name.test<-c('p.fisher.less','prior')

      ans.tmp<-rep(NA,length(name.test))
      sum1.A<-sum(x1[1:n1])
      if(sum1.A<=n.stop1)
      {
        ans.A<-'A.fail.stage1'
        p1.est<-mean(x1)
      }else
      {
        sum2.A<-sum(x1)
        if(sum2.A<=n.stop2)
        {
          ans.A<-'A.fail.stage2'
          p1.est<-mean(x1)
        } else
        {
          ans.A<-'A.pass'
          p1.est<-mean(x1)
        }
      }

      sum1.B<-sum(x2[1:n1])
      if(sum1.B<=n.stop1)
      {
        ans.B<-'B.fail.stage1'
        p2.est<-mean(x2)
      } else
      {
        sum2.B<-sum(x2)
        if(sum2.B<=n.stop2)
        {
          ans.B<-'B.fail.stage2'
          p2.est<-mean(x2)
        } else
        {
          ans.B<-'B.pass'
          p2.est<-mean(x2)
        }
      }

      if(ans.A=='A.pass' & ans.B=='B.pass')
      {
        table1<-cbind(table(factor(x1,level=c(1,0))),table(factor(x2,level=c(1,0))))

        ans.tmp<-fun.bay.fisher.v(as.vector(table1))
        names(ans.tmp)<-name.test
      }
      ans.tmp1<-c(ans.A,ans.B,sum(x1),length(x1)-sum(x1),sum(x2),length(x2)-sum(x2),ans.tmp)
      names(ans.tmp1)<-c('arm.A','arm.B','R_A','NR_A','R_B','NR_B',name.test)
      ans<-rbind(ans,ans.tmp1)
    }
    ans
  }

  baye.sim.single.stage.fun<-function(sim.n0,p1.true,p2.true, n.total,n.stop)
  {
    ans<-numeric(0)

    for(i in 1:sim.n0)
    {
      x1<-rbinom(n.total,1,p1.true)
      x2<-rbinom(n.total,1,p2.true)

      ans.A<-ans.B<-NULL
      p.tmp<-poster.mean1<-NA
      name.test<-c('p.fisher.less','prior')

      ans.tmp<-rep(NA,length(name.test))
      sum.A<-sum(x1)
      if(sum.A<=n.stop)
      {
        ans.A<-'A.fail'
        p1.est<-mean(x1)
      }else
        {
          ans.A<-'A.pass'
          p1.est<-mean(x1)
        }


      sum.B<-sum(x2)
      if(sum.B<=n.stop)
      {
        ans.B<-'B.fail'
        p2.est<-mean(x2)
      } else
        {
          ans.B<-'B.pass'
          p2.est<-mean(x2)
        }

      if(ans.A=='A.pass' & ans.B=='B.pass')
      {
        table1<-cbind(table(factor(x1,level=c(1,0))),table(factor(x2,level=c(1,0))))

        ans.tmp<-fun.bay.fisher.v(as.vector(table1))
        names(ans.tmp)<-name.test
      }
      ans.tmp1<-c(ans.A,ans.B,sum(x1),length(x1)-sum(x1),sum(x2),length(x2)-sum(x2),ans.tmp)
      names(ans.tmp1)<-c('arm.A','arm.B','R_A','NR_A','R_B','NR_B',name.test)
      ans<-rbind(ans,ans.tmp1)
    }
    ans
  }

  p0.true.list<-seq(p10,p20,by=.05)
  if(length(p0.true.list)>3)
  {
    p1.true.list<-seq(p10,p20,by=.05)[1:3]
    p2.true.list<-rev(seq(p20,p10,by=-0.05)[1:3])
  } else {
    p1.true.list<-seq(p10,p20,len=5)[1:3]
    p2.true.list<-rev(seq(p20,p10,len=5)[1:3])

  }
  my.p1p2.list<-list(p10=p10,p20=p20, p1.true.list=p1.true.list,p2.true.list=p2.true.list)


  sim.ans<-ph2simon(p10, p20, alpha1, beta1)
  if(optimal) tmp1<-sim.ans$out[,'EN(p0)'] else tmp1<-sim.ans$out[,'n']

  nn<-sim.ans$out[tmp1==min(tmp1),1:4]
  n.total<-nn['n']
  n1<-nn['n1']
  n.stop1<-nn['r1']
  n.stop2<-nn['r']
  ans.list<-list()

  kk0<-1
  ans<-baye.sim.fun(sim.n0=sim.n0,n.total=n.total,n1=n1,p1.true=p10,p2.true=p10,n.stop1=n.stop1,n.stop2=n.stop2)
  ans.list[[kk0]]<-ans
  names(ans.list)[kk0]<-paste('p1_',p10,'_p2_',p10,sep='_')

  p1.p2.table<-c(p10,p10)

  for(kk1 in 1:length(p1.true.list))
  {
    p1.true<-p1.true.list[kk1]
    for(kk2 in kk1:length(p2.true.list))
    {
      kk0<-kk0+1
      p2.true<-p2.true.list[kk2]
      ans<-baye.sim.fun(sim.n0=sim.n0,n.total=n.total,n1=n1,p1.true=p1.true,p2.true=p2.true,n.stop1=n.stop1,n.stop2=n.stop2)

      ans.list[[kk0]]<-ans
      names(ans.list)[kk0]<-paste('p1_',p1.true,'_p2_',p2.true,sep='_')
      p1.p2.table<-rbind(p1.p2.table,c(p1.true,p2.true))

    } #    for(kk1 in 1:length(p1.true.list))

  } # for(kk2 in 1:length(p2.true.list))

  n.design<-sim.ans$out[tmp1==min(tmp1),,drop=F]
  dimnames(n.design)[[1]][1]<-study.design
  my.p1p2.list$design<-study.design
  my.p1p2.list$n.max<-n.max
  dimnames(p1.p2.table)[[1]]<-names(ans.list)
  list(ans.list=ans.list,sim.ans=sim.ans,my.p1p2.list=my.p1p2.list,n.design=n.design,p1.p2.table=p1.p2.table)
}


Bayesian.simon.single.stage.fun<-function(sim.n0=10,  beta1=.1,alpha1=.1, p10=0.2,p20=0.4,n.max=NULL)
{
  require(clinfun)
#  set.seed(1000)

  prob.mean.v<-function(data)
  {
    sim.n<-100000
    aa1<-data[1] # number of response in Arm A
    bb1<-data[2] # number of non-response in Arm A
    aa2<-data[3] # number of response in Arm B
    bb2<-data[4] # number of non-response in Arm B
    beta1<-rbeta(sim.n,aa1,bb1);
    beta2<-rbeta(sim.n,aa2,bb2)
    mean(beta1<beta2)
  }

  fun.bay.fisher.v<-function(x,prior_armA=c(1,1),prior_armB=c(1,1)){
    table1<-cbind(x[1:2],x[3:4])
    # x[1:2]= (# of response, # of non-response) in arm A
    # x[3:4]= (# of response, # of non-response) in arm B
    # prior_armA= two beta parameters: a and b for Arm A
    # prior_armB= two beta parameters: a and b for Arm B

    p.fisher.less<-fisher.test(table1,alternative ='less')$p.value
    c(p.fisher.less,prob.mean.v(x+c(prior_armA,prior_armB)))
  }


  baye.sim.single.stage.fun<-function(sim.n0,p1.true,p2.true, n.total,n.stop)
  {
    ans<-numeric(0)

    for(i in 1:sim.n0)
    {
      x1<-rbinom(n.total,1,p1.true)
      x2<-rbinom(n.total,1,p2.true)

      ans.A<-ans.B<-NULL
      p.tmp<-poster.mean1<-NA
      name.test<-c('p.fisher.less','prior')

      ans.tmp<-rep(NA,length(name.test))
      sum.A<-sum(x1)
      if(sum.A<=n.stop)
      {
        ans.A<-'A.fail'
        p1.est<-mean(x1)
      }else
      {
        ans.A<-'A.pass'
        p1.est<-mean(x1)
      }


      sum.B<-sum(x2)
      if(sum.B<=n.stop)
      {
        ans.B<-'B.fail'
        p2.est<-mean(x2)
      } else
      {
        ans.B<-'B.pass'
        p2.est<-mean(x2)
      }

      if(ans.A=='A.pass' & ans.B=='B.pass')
      {
        table1<-cbind(table(factor(x1,level=c(1,0))),table(factor(x2,level=c(1,0))))

        ans.tmp<-fun.bay.fisher.v(as.vector(table1))
        names(ans.tmp)<-name.test
      }
      ans.tmp1<-c(ans.A,ans.B,sum(x1),length(x1)-sum(x1),sum(x2),length(x2)-sum(x2),ans.tmp)
      names(ans.tmp1)<-c('arm.A','arm.B','R_A','NR_A','R_B','NR_B',name.test)
      ans<-rbind(ans,ans.tmp1)
    }
    ans
  }

  p0.true.list<-seq(p10,p20,by=.05)
  if(length(p0.true.list)>3)
  {
    p1.true.list<-seq(p10,p20,by=.05)[1:3]
    p2.true.list<-rev(seq(p20,p10,by=-0.05)[1:3])
  } else {
    p1.true.list<-seq(p10,p20,len=5)[1:3]
    p2.true.list<-rev(seq(p20,p10,len=5)[1:3])

  }
  my.p1p2.list<-list(p10=p10,p20=p20, p1.true.list=p1.true.list,p2.true.list=p2.true.list)


  sim.ans<-ph2single(p10, p20, alpha1, beta1)
  tmp1<-sim.ans[1,]

  n.total<-tmp1[1,'n']
  n.stop<-tmp1[1,'r']
  ans.list<-list()

  kk0<-1
  ans<-baye.sim.single.stage.fun(sim.n0=sim.n0,n.total=n.total,p1.true=p10,p2.true=p10,n.stop=n.stop)
  ans.list[[kk0]]<-ans
  names(ans.list)[kk0]<-paste('p1_',p10,'_p2_',p10,sep='_')

  p1.p2.table<-c(p10,p10)

  for(kk1 in 1:length(p1.true.list))
  {
    p1.true<-p1.true.list[kk1]
    for(kk2 in kk1:length(p2.true.list))
    {
      kk0<-kk0+1
      p2.true<-p2.true.list[kk2]
      ans<-baye.sim.single.stage.fun(sim.n0=sim.n0,n.total=n.total,p1.true=p1.true,p2.true=p2.true,n.stop=n.stop)

      ans.list[[kk0]]<-ans
      names(ans.list)[kk0]<-paste('p1_',p1.true,'_p2_',p2.true,sep='_')
      p1.p2.table<-rbind(p1.p2.table,c(p1.true,p2.true))

    } #    for(kk1 in 1:length(p1.true.list))

  } # for(kk2 in 1:length(p2.true.list))

  n.design<-sim.ans[1,]
  dimnames(n.design)[[1]][1]<-'single stage design'
  my.p1p2.list$design<-'single stage design'
  my.p1p2.list$n.max<-n.max
  dimnames(p1.p2.table)[[1]]<-names(ans.list)
  list(ans.list=ans.list,sim.ans=sim.ans,my.p1p2.list=my.p1p2.list,n.design=n.design,p1.p2.table=p1.p2.table)
}




ui <- fluidPage(
  headerPanel('A Bayesian Pick-the-Winner Design in a Randomized Phase II Clinical Trial'),

  tabsetPanel(

    tabPanel("Baysian pick-the-winner design (two-stage)",
             sidebarLayout(
               sidebarPanel(

                 sliderInput("alpha",label="Type I error",min = 0, max = 1, value = 0.1),

                 sliderInput("beta",label="Type II error",min = 0, max = 1, value = 0.1),

                 sliderInput("p10",label="Response rate in Arm A",min = 0, max = 1, value =0.2),

                 sliderInput("p20",label="Response rate in Arm B",min = 0, max = 1, value =0.4),

                 selectInput("design",label="Study design",choices=c('optimal','MiniMax')),

                 numericInput("num", label ="Number of Simulations", value = 100),

                 actionButton("Submit2","Calculate"),
                 downloadButton('downloadReport')

               ),
               mainPanel(

#                 h2("An integrated Bayesian posterior probability with Simon two-stage design for a randomized phase II clinical trial"),
#                 tags$hr(),

                 h4("Hypothesis for Power Analysis"),
                 textOutput("hypothesis"),
                 tags$hr(),

                 h4("Sample Size Calculation"),
                 textOutput("samplesize"),
                 tags$hr(),

                 h4("Operating Characteristics"),
                 textOutput("operation"),
                 tags$hr(),

                 h4("Power Analysis"),
                 textOutput("power"),
                 tags$hr(),

                 h4("Type I error"),
                 textOutput("typeI"),
                 tags$hr(),

                 h4("Summary"),
                 textOutput("summary"),
                 tags$hr(),

                 h4("Table of Power Analysis"),
                 verbatimTextOutput("table")



                 ))),

tabPanel("Baysian pick-the-winner design (single stage)",
         sidebarLayout(
           sidebarPanel(

             sliderInput("alpha",label="Type I error",min = 0, max = 1, value = 0.1),

             sliderInput("beta",label="Type II error",min = 0, max = 1, value = 0.1),

             sliderInput("p10",label="Response rate in Arm A",min = 0, max = 1, value =0.2),

             sliderInput("p20",label="Response rate in Arm B",min = 0, max = 1, value =0.4),

             numericInput("num", label ="Number of Simulations", value = 100),

             actionButton("SubmitSingleStage","Calculate")

           ),
           mainPanel(

             h4("Hypothesis for Power Analysis"),
             textOutput("hypothesis1"),
             tags$hr(),

             h4("Sample Size Calculation"),
             textOutput("samplesize1"),
             tags$hr(),

             h4("Operating Characteristics"),
             textOutput("operation1"),
             tags$hr(),

             h4("Power Analysis"),
             textOutput("power1"),
             tags$hr(),

             h4("Type I error"),
             textOutput("typeI1"),
             tags$hr(),

             h4("Summary"),
             textOutput("summary1"),
             tags$hr(),

             h4("Table of Power Analysis"),
             verbatimTextOutput("table1")

           ))),

tabPanel("Calculaiton of Bayesian posterior probability",
         sidebarLayout(
           sidebarPanel(

             numericInput("numArmAResponse", label ="Number of response in arm A", value = 10),
             numericInput("numArmANoResponse", label ="Number of non-response in arm A", value = 10),
             numericInput("numArmBResponse", label ="Number of response in arm B", value = 10),
             numericInput("numArmBNoResponse", label ="Number of non-response in arm B", value = 10),
             actionButton("SubmitBayeProb","Calculate")
           ),
           mainPanel(

             h4("Bayesian posterior probability"),
             plotOutput("BayeProbPlot")

           )))
  )

  )


server <- function(input,output){

  #---two stage stage--------

  get.result <- eventReactive(input$Submit2, {

    sim.n0.tmp=input$num
    study.design.tmp=input$design
    alpha1.tmp=input$alpha
    beta1.tmp=input$beta
    p10.tmp=input$p10
    p20.tmp=input$p20

    tmp99<-Bayesian.simon.two.stage.fun(
      sim.n0=sim.n0.tmp, study.design=study.design.tmp,beta1=beta1.tmp,alpha1=alpha1.tmp,
      p10=p10.tmp,p20=p20.tmp)

    p10<-tmp99$my.p1p2.list$p10
    p20<-tmp99$my.p1p2.list$p20
    alpha1<-tmp99$sim.ans$alpha
    beta1<-tmp99$sim.ans$beta
    n.max<-tmp99$my.p1p2.list$n.max
    ans.list.tmp<-tmp99$ans.list
    design.status<-tmp99$my.p1p2.list$design
    p1.p2.table<-tmp99$p1.p2.table
    p1.p2.diff<-round(apply(p1.p2.table,1,diff),10)
    p1.p2.diff.uni<-sort(unique(p1.p2.diff))
    scenario<-names(tmp99$ans.list)
    if(design.status=='optimal')  design.name<-'Optimal' else design.name<-'Mini-Max'

    nn<-tmp99$n.design
    n.total<-nn[1,'n']
    n1<-nn[1,'n1']
    n.stop1<-nn[1,'r1']
    n.stop2<-nn[1,'r']
    ans.list<-ans.list.tmp[scenario]
    bb1<-sapply(ans.list,function(ans) {
      sim.n0<-dim(ans)[1]
      ans.A<-factor(ans[,'arm.A'],level=c("A.fail.stage1", "A.fail.stage2", "A.pass"))
      ans.B<-factor(ans[,'arm.B'],level=c("B.fail.stage1", "B.fail.stage2", "B.pass"))
      table1<-table(ans.A,ans.B)/sim.n0
      bay.prob<-as.numeric(ans[,'prior'])
      winner.prob<-sum(bay.prob>0.8,na.rm=T)/sim.n0
      list(table=round(table1,2),
           Bay.prob.B.winner=winner.prob,
           prob.B.pass.A.fail=sum(table1[1:2,3]),
           power.B=sum(table1[1:2,3])+winner.prob)
    }
    ,simplify=F)



    a1<-paste('From historical data, we will consider ',round(p10*100),'% response rate as not warranting further study. We will use ',round(p20*100),'% response rate as a promising result to pursue further study. In other words, we are interested in at least ',round((p20-p10)*100),'% (',round(p20*100),'% vs. ',round(p10*100),'%) improvement in treatment efficacy for arms B versus A. For each arm, using a Simon ',design.name, ' two-stage design with ',round(alpha1*100),'% type I error rate and ', round(beta1*100),'% type II error rate, ',n1 ,' patients will be enrolled in the first stage of the trial. If ',n.stop1 ,' or fewer patients respond, the treatment will be stopped. If ',n.stop1+1 ,' or more patients show a response, ', n.total-n1,' additional patients (a total of ',n.total ,' patients per group) will be enrolled. If the total number responding is ',n.stop2,' or less, we will conclude that the treatment is not effective.',if(!is.null(n.max)) paste(' We plan to enroll a maximum of ',n.max, ' patients in order to have ',round( (n.total*2)/n.max*100), '% of the enrolled patients remain on the trial if both arms finish the 2nd stage.',sep=''),' If both arms fail at the first or second stage, the trial will stop. No winner will be claimed. The sample size will be ',n1*2, ' if both arms fail at the first stage and ',n1+n.total, ' if only one arm fails at the first stage. If only one arm pass the second stage, the arm will be the winner. If both arms pass the second stage, we will use the posterior probability, $Pr(B>A)$, (probability of the response rate in arm B higher than in arm A) to select the winner. A non-informative prior of beta distribution, beta(1,1) in both arms will be used to calculate the posterior probability. Arm B will be claimed as the winner if $Pr(B>A)>\\delta=$ 0.8.',sep='')

    a2<-paste('The operating characteristics of the design is evaluated by simulation (',sim.n0.tmp,' times) using R software (www.r-project.org) with "clinfun" package. In particular, we are interested in the probability of (correctly) selecting an arm as superior to the other arm if it is truly superior, and conversely, the probability of (incorrectly) selecting an arm that is no better than the other arm. ',sep='')

    bb1.largest.diff<-bb1[p1.p2.diff==p1.p2.diff.uni[4]][[1]]
    bb1.2nd_largest.diff<-bb1[p1.p2.diff==p1.p2.diff.uni[3]]
    p1.p2.table_2nd_largest.diff<-p1.p2.table[p1.p2.diff==p1.p2.diff.uni[3],]
    bb1.3rd_largest.diff<-bb1[p1.p2.diff==p1.p2.diff.uni[2]]
    p1.p2.table_3rd_largest.diff<-p1.p2.table[p1.p2.diff==p1.p2.diff.uni[2],]

    bb1.no.diff<-bb1[p1.p2.diff==p1.p2.diff.uni[1]][[1]]

    a21<-paste('Power: Assuming that the true probabilities of response in arms B and A are ',p20*100,'% and ',p10*100,'%, respectively (scenario 1: ',round((p20-p10)*100),'% difference of response rate), the overall probability (power) of correctly choosing arm B as superior is ',round(bb1.largest.diff$power.B*100),'% on the basis of superiority shown at the end of the trial. The probability of stopping arm A early and declaring arm B superior at the end of the trial is ', round(sum(bb1.largest.diff$prob.B.pass.A.fail)*100),'%. There are ',round(sum(bb1.largest.diff$table['A.pass','B.pass'])*100),'% of both arms passing the second stage with ',round(bb1.largest.diff$Bay.prob.B.winner*100),'% claiming arm B as the winner by the Bayesian posterior probability. In a ',round(p1.p2.diff.uni[3]*100),'% difference of response rate, the overall power is ',round(bb1.2nd_largest.diff[[1]]$power.B*100),'% and ',round(bb1.2nd_largest.diff[[2]]$power.B*100),'% for the comparison of arms B and A with ',round(p1.p2.table_2nd_largest.diff[1,2]*100),'% versus ',round(p1.p2.table_2nd_largest.diff[1,1]*100),'% (scenario 2) and ',round(p1.p2.table_2nd_largest.diff[2,2]*100),'% versus ',round(p1.p2.table_2nd_largest.diff[2,1]*100),'% (scenario 3), respectively. Proportion of both arms passing the 2nd stage is ',round(bb1.2nd_largest.diff[[1]]$table['A.pass','B.pass']*100),'% in scenario 2 (scenario 3: ',round(bb1.2nd_largest.diff[[2]]$table['A.pass','B.pass']*100),'%), with ',round(bb1.2nd_largest.diff[[1]]$Bay.prob.B.winner*100),'% (scenario 3: ',round(bb1.2nd_largest.diff[[2]]$Bay.prob.B.winner*100),'%) claiming arm B as the winner by the Bayesian posterior probability. ',sep='')

    a22<-paste('Type I error: In the null hypothesis of a ',p10*100,'% response rate in both arms, there are ',round(bb1.no.diff$power.B*100),'% misclassifying arm B as winner (i.e., ',round(bb1.no.diff$power.B*100),'% type I error). Among them, only ',round(bb1.no.diff$table['A.pass','B.pass']*100,2),'% has both arms passing the 2nd stage, and less than ',round(bb1.no.diff$Bay.prob.B.winner*100,2),'% misclassify arm B as winner. ',sep='')
    a23<-paste('Summary: With $\\delta$=0.8, the design has a ',round(bb1.largest.diff$power.B*100),'% power to detect a ',round((p20-p10)*100),'% difference of response rate. The power decreases to a range of ',paste(range(sapply(bb1.2nd_largest.diff,function(x) round(x$power.B*100))),collapse='-'),'% to differentiate a ',round(p1.p2.diff.uni[3]*100),'% difference of response rate. The type I error is controlled at ',round(bb1.no.diff$power.B*100),'% when both arms have a ',p10*100,'% response rate. ',sep='')


    fun1<-function(x)
    {
      rbind(x$table,c(paste('Both arms passing the 2nd stage: ', round(x$table['A.pass','B.pass']*100,2),'%. Among them,  Arm B claims ',round(x$Bay.prob.B.winner*100,2),'% as winner',sep=''),paste('Overall power of Arm B= ',round(x$power.B*100),'%',sep=''), paste('Type I error= ', round(x$power.B*100,2),'%',sep='')))
    }
    fun2<-function(x)
    {
      tmp0<-sub('p2__','Arm B=',sub('p1__','Arm A=',x))
      paste('Scenario :',sub('.*__','',tmp0),' versus ',sub('__.*','',tmp0),sep='')
    }

    bb1.selected<-c(bb1[p1.p2.diff==p1.p2.diff.uni[4]],bb1[p1.p2.diff==p1.p2.diff.uni[3]],bb1[p1.p2.diff==p1.p2.diff.uni[1]])


    tmp10<-sapply(bb1.selected,fun1,simplify=F)
    tmp20<-numeric()
    name1<-names(tmp10)
    tmp30<-numeric()
    for(i in 1:length(name1))
    {
      tmp30<-c(tmp30,c(sub('Scenario :',paste('Scenario ',i, ': ',sep=''), fun2(name1[i])),rep('',3)))
    }
    for(i in 1:length(tmp10))
      tmp20<-rbind(tmp20,tmp10[[i]])
    tmp20<-cbind(tmp30,tmp20)
    tmp98<-list(ans=tmp20,tmp=tmp10,a1=a1,a2=a2,power=a21,typeI=a22,a3=a23,p10=p10,p20=p20)

    return(tmp98)


  })


  output$hypothesis <- renderText({
    res <- get.result()

    paste('Comparison of ',round(res$p20*100),'% versus ',round(res$p10*100),'% Response Rate')

  })

  output$samplesize <- renderText({
    get.result()$a1
  })

  output$operation <- renderText({
    get.result()$a2
  })

  output$power <- renderText({
    get.result()$power
  })

  output$typeI <- renderText({
    get.result()$typeI
  })

  output$summary <- renderText({
    get.result()$a3
  })

  output$table <- renderPrint({

    tmp98 <- get.result()

    tmp2<-as.vector(tmp98$ans[,1])
    name1<-tmp2[tmp2!='']
    #name1<-paste('Scenario ',1:length(name1),sub('Scenario :',': ',name1),sep='')
    tmp3<-tmp98$tmp

    for(i in 1:length(tmp3))
    {
      tmp40<-tmp3[[i]]
      tmp4<-data.frame(tmp40[1:3,])
      tmp41<-as.vector(tmp40[4,])
      name2<- if(i!=4) paste(name1[i],' (',tmp41[2],')',sep='') else paste(name1[i],' (',tmp41[3],')',sep='')

      cat('\n-----------------------------------------------\n')
      cat('\n\n',name2,'\n')
      print(kable(tmp4))
      cat('\n',tmp41[1],'\n')
      if(i!=4) cat('\n',tmp41[2],'\n') else cat('\n',tmp41[3],'\n')
    }

  })

  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('BayesianPickWinnerReport.docx')
    },

    content = function(file) {
      out <- render('BayesianPickWinnerReport.Rmd')
      file.rename(out, file)
    }
  )
  #---single stage--------

  get.result.single.stage <- eventReactive(input$SubmitSingleStage, {

    sim.n0.tmp=input$num
    study.design.tmp=input$design
    alpha1.tmp=input$alpha
    beta1.tmp=input$beta
    p10.tmp=input$p10
    p20.tmp=input$p20

    tmp99<-Bayesian.simon.single.stage.fun(
      sim.n0=sim.n0.tmp, beta1=beta1.tmp,alpha1=alpha1.tmp,
      p10=p10.tmp,p20=p20.tmp)


    p10<-tmp99$my.p1p2.list$p10
    p20<-tmp99$my.p1p2.list$p20
    alpha1<-alpha1.tmp
    beta1<-beta1.tmp
    n.max<-tmp99$my.p1p2.list$n.max
    ans.list.tmp<-tmp99$ans.list
    design.status<-tmp99$my.p1p2.list$design
    p1.p2.table<-tmp99$p1.p2.table
    p1.p2.diff<-round(apply(p1.p2.table,1,diff),10)
    p1.p2.diff.uni<-sort(unique(p1.p2.diff))
    scenario<-names(tmp99$ans.list)
    design.name<-tmp99$my.p1p2.list$design

    nn<-tmp99$n.design
    n.total<-nn[1,'n']
    n.stop<-nn[1,'r']
    ans.list<-ans.list.tmp[scenario]
    bb1<-sapply(ans.list,function(ans) {
      sim.n0<-dim(ans)[1]
      ans.A<-factor(ans[,'arm.A'],level=c("A.fail",  "A.pass"))
      ans.B<-factor(ans[,'arm.B'],level=c("B.fail",  "B.pass"))
      table1<-table(ans.A,ans.B)/sim.n0
      bay.prob<-as.numeric(ans[,'prior'])
      winner.prob<-sum(bay.prob>0.8,na.rm=T)/sim.n0
      list(table=round(table1,2),
           Bay.prob.B.winner=winner.prob,
           prob.B.pass.A.fail=sum(table1[1,2]),
           power.B=sum(table1[1,2])+winner.prob)
    }
    ,simplify=F)

    a1<-paste('From historical data, we will consider ',round(p10*100),'% response rate as not warranting further study. We will use ',round(p20*100),'% response rate as a promising result to pursue further study. In other words, we are interested in at least ',round((p20-p10)*100),'% (',round(p20*100),'% vs. ',round(p10*100),'%) improvement in treatment efficacy for arms B versus A. For each arm, using the Fleming single stage design with ',round(alpha1*100),'% type I error rate and ', round(beta1*100),'% type II error rate, ',n.total ,' patients will be enrolled in the trial. If ',n.stop+1 ,' patients or more respond, the treatment will be considered competitive. ', if(!is.null(n.max)) paste(' We plan to enroll a maximum of ',n.max, ' patients in order to have ',round( (n.total*2)/n.max*100), '% of the enrolled patients remain on the trial if both arms finish at end of the trial.',sep=''),' If both arms fail, the trial will stop. No winner will be claimed. The sample size will be ',n.total, ' patients per arm. If only one arm is competitive  (i.e., number of responses >',n.stop,'), the arm will be the winner. If both arms are competitive, we will use the posterior probability, $Pr(B>A)$, (probability of the response rate in arm B higher than in arm A) to select the winner. A non-informative prior of beta distribution, beta(1,1), in both arms will be used to calculate the posterior probability. Arm B will be claimed as the winner if $Pr(B>A)>\\delta=$ 0.8.',sep='')

    a2<-paste('The operating characteristics of the design is evaluated by simulation (',sim.n0.tmp,' times) using R software (www.r-project.org) with "clinfun" package. In particular, we are interested in the probability of (correctly) selecting an arm as superior to the other arm if it is truly superior, and conversely, the probability of (incorrectly) selecting an arm that is no better than the other arm. ',sep='')

    bb1.largest.diff<-bb1[p1.p2.diff==p1.p2.diff.uni[4]][[1]]
    bb1.2nd_largest.diff<-bb1[p1.p2.diff==p1.p2.diff.uni[3]]
    p1.p2.table_2nd_largest.diff<-p1.p2.table[p1.p2.diff==p1.p2.diff.uni[3],]
    bb1.3rd_largest.diff<-bb1[p1.p2.diff==p1.p2.diff.uni[2]]
    p1.p2.table_3rd_largest.diff<-p1.p2.table[p1.p2.diff==p1.p2.diff.uni[2],]

    bb1.no.diff<-bb1[p1.p2.diff==p1.p2.diff.uni[1]][[1]]

    a21<-paste('Power: Assuming that the true probabilities of response in arms B and A are ',p20*100,'% and ',p10*100,'%, respectively (scenario 1: ',round((p20-p10)*100),'% difference of response rate), the overall probability (power) of correctly choosing arm B as superior is ',round(bb1.largest.diff$power.B*100),'% on the basis of superiority shown at the end of the trial. The probability of stopping arm A early and declaring arm B superior at the end of the trial is ', round(sum(bb1.largest.diff$prob.B.pass.A.fail)*100),'%. There are ',round(sum(bb1.largest.diff$table['A.pass','B.pass'])*100),'% of both arms being competitive with ',round(bb1.largest.diff$Bay.prob.B.winner*100),'% claiming arm B as the winner by the Bayesian posterior probability. In a ',round(p1.p2.diff.uni[3]*100),'% difference of response rate, the overall power is ',round(bb1.2nd_largest.diff[[1]]$power.B*100),'% and ',round(bb1.2nd_largest.diff[[2]]$power.B*100),'% for the comparison of arms B and A with ',round(p1.p2.table_2nd_largest.diff[1,2]*100),'% versus ',round(p1.p2.table_2nd_largest.diff[1,1]*100),'% (scenario 2) and ',round(p1.p2.table_2nd_largest.diff[2,2]*100),'% versus ',round(p1.p2.table_2nd_largest.diff[2,1]*100),'% (scenario 3), respectively. Proportion of both arms being competitive  is ',round(bb1.2nd_largest.diff[[1]]$table['A.pass','B.pass']*100),'% in scenario 2 (scenario 3: ',round(bb1.2nd_largest.diff[[2]]$table['A.pass','B.pass']*100),'%), with ',round(bb1.2nd_largest.diff[[1]]$Bay.prob.B.winner*100),'% (scenario 3: ',round(bb1.2nd_largest.diff[[2]]$Bay.prob.B.winner*100),'%) claiming arm B as the winner by the Bayesian posterior probability. ',sep='')

    a22<-paste('Type I error: In the null hypothesis of a ',p10*100,'% response rate in both arms, there are ',round(bb1.no.diff$power.B*100),'% misclassifying arm B as winner (i.e., ',round(bb1.no.diff$power.B*100),'% type I error). Among them, only ',round(bb1.no.diff$table['A.pass','B.pass']*100,2),'% has both arms being competitive, and ',round(bb1.no.diff$Bay.prob.B.winner*100,2),'% misclassify arm B as winner. ',sep='')
    a23<-paste('Summary: With $\\delta$=0.8, the design has a ',round(bb1.largest.diff$power.B*100),'% power to detect a ',round((p20-p10)*100),'% difference of response rate. The power decreases to a range of ',paste(range(sapply(bb1.2nd_largest.diff,function(x) round(x$power.B*100))),collapse='-'),'% to differentiate a ',round(p1.p2.diff.uni[3]*100),'% difference of response rate. The type I error is controlled at ',round(bb1.no.diff$power.B*100),'% when both arms have a ',p10*100,'% response rate. ',sep='')


    fun1<-function(x)
    {
      rbind(x$table,c(paste('Both arms pass: ', round(x$table['A.pass','B.pass']*100,2),'%. Among them,  Arm B claims ',round(x$Bay.prob.B.winner*100,2),'% as winner',sep=''),paste('Overall power of Arm B= ',round(x$power.B*100),'%',sep='')))
    }
    fun2<-function(x)
    {
      tmp0<-sub('p2__','Arm B=',sub('p1__','Arm A=',x))
      paste('Scenario :',sub('.*__','',tmp0),' versus ',sub('__.*','',tmp0),sep='')
    }

    bb1.selected<-c(bb1[p1.p2.diff==p1.p2.diff.uni[4]],bb1[p1.p2.diff==p1.p2.diff.uni[3]],bb1[p1.p2.diff==p1.p2.diff.uni[1]])


    tmp10<-sapply(bb1.selected,fun1,simplify=F)
    tmp20<-numeric()
    name1<-names(tmp10)
    tmp30<-numeric()
    for(i in 1:length(name1))
    {
      tmp30<-c(tmp30,c(sub('Scenario :',paste('Scenario ',i, ': ',sep=''), fun2(name1[i])),rep('',2)))
    }
    for(i in 1:length(tmp10))
      tmp20<-rbind(tmp20,tmp10[[i]])
    tmp20<-cbind(tmp30,tmp20)
    tmp98<-list(ans=tmp20,tmp=tmp10,a1=a1,a2=a2,power=a21,typeI=a22,a3=a23,p10=p10,p20=p20)

    return(tmp98)


  })


  output$hypothesis1 <- renderText({
    res <- get.result.single.stage()

    paste('Comparison of ',round(res$p20*100),'% versus ',round(res$p10*100),'% Response Rate')

  })

  output$samplesize1 <- renderText({
    get.result.single.stage()$a1
  })

  output$operation1 <- renderText({
    get.result.single.stage()$a2
  })

  output$power1 <- renderText({
    get.result.single.stage()$power
  })

  output$typeI1 <- renderText({
    get.result.single.stage()$typeI
  })

  output$summary1 <- renderText({
    get.result.single.stage()$a3
  })

  output$table1 <- renderPrint({

    tmp98 <- get.result.single.stage()

    tmp2<-as.vector(tmp98$ans[,1])
    name1<-tmp2[tmp2!='']
    #name1<-paste('Scenario ',1:length(name1),sub('Scenario :',': ',name1),sep='')
    tmp3<-tmp98$tmp

    for(i in 1:length(tmp3))
    {
      tmp40<-tmp3[[i]]
      tmp4<-data.frame(tmp40[1:2,])
      tmp41<-as.vector(tmp40[3,])
      name2<- if(i!=4) paste(name1[i],' (',tmp41[2],')',sep='') else paste(name1[i],' (',sub('Overall power of Arm B','Type I error',tmp41[2]),')',sep='')

      cat('\n-----------------------------------------------\n')
      cat('\n\n',name2,'\n')
      print(kable(tmp4))
      cat('\n',tmp41[1],'\n')
      if(i!=4) cat('\n',tmp41[2],'\n') else cat('\n',sub('Overall power of Arm B','Type I error',tmp41[2]),'\n')
    }

  })



  #---this part is for Bayesian posterior prob.

  get.result.baye.prob <- eventReactive(input$SubmitBayeProb, {

    tmp98<-Bayesian_posterior_probability(n_response_armA=input$numArmAResponse, n_nonresponse_armA=input$numArmANoResponse, n_response_armB=input$numArmBResponse, n_nonresponse_armB=input$numArmBNoResponse,sim.n=100000)
    list(prob=tmp98,n_response_armA=input$numArmAResponse, n_nonresponse_armA=input$numArmANoResponse, n_response_armB=input$numArmBResponse, n_nonresponse_armB=input$numArmBNoResponse)
  })


  output$BayeProbPlot <- renderPlot({
    seq1<-seq(0,1,len=10000)
    res<-get.result.baye.prob()
    armA.a<-1 + res$n_response_armA
    armA.b<-1 + res$n_nonresponse_armA
    armB.a<-1 + res$n_response_armB
    armB.b<-1 + res$n_nonresponse_armB
    density.armA<-dbeta(seq1,armA.a,armA.b)
    density.armB<-dbeta(seq1,armB.a,armB.b)
    plot(range(seq1),range(c(density.armA,density.armB)),xlab='response rate',ylab='density',type='n')
    legend(min(seq1),max(c(density.armA,density.armB)),c('arm A','arm B'),col=1:2,lty=2:1)
    lines(seq1,density.armA,lty=2)
    lines(seq1,density.armB,col=2,lty=1)
    title(paste('Bayesian Posterior Probability (Pr(arm B>arm A))=',round(res$prob,4),'\n (non-informative beta prior, beta(1,1), in each arm)',sep=''))
  })

}

shinyApp(ui=ui,server=server)
