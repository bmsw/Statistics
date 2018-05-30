

flowDiv <- function (myworkspaces=list(), gate.name=NULL, channel1=NULL, channel2=NULL, nbins=16, dil=c(), flowrate=c(), volume=c(), nsamples=NULL, use.beads=FALSE, beads.gate=NULL, doplot=TRUE, file.name=NULL, save.csv=FALSE, ialpha="invsimpson", ibeta="bray", runalpha="TRUE", runbeta="TRUE"){


  samples=list()
  MatList = list()
  id.files=c()
  mean_chan1=c()
  mean_chan2=c()
  timediff.all=c()
  total.flow = c()
  timediff.all=c()
  if(is.null(dil)) dil=rep(1, nsamples)

  for(i in 1:length(myworkspaces)){




    workspace<-flowWorkspace::openWorkspace(myworkspaces[[i]])


    gating.set<-flowWorkspace::parseWorkspace(workspace, name=1)

    flowList = list()
    matrixList = list()
    id = c()
    beadsList = list()
    fmean_chan1=c()
    fmean_chan2=c()
    c1.invTrans = list()
    c2.invTrans = list()
    versions=c()
    invstate="off"


    etim=c()
    btim=c()
    timediff=c()


    for (z in 1:length(gating.set))
    {

      id[z] = sub("\\.fcs", "", tail(unlist(strsplit(flowWorkspace::getData(gating.set[[z]])@description$FILENAME, split="\\/")), n=1))



      if(!is.null(flowWorkspace::getTransformations(gating.set[[z]], channel=channel1, inverse=TRUE)))
      {
        c1.invTrans[[z]] = flowWorkspace::getTransformations(gating.set[[z]], channel=channel1, inverse=TRUE)
        invstate="on"
      }
      else  c1.invTrans[[z]] = flowWorkspace::getTransformations(gating.set[[z]], channel=channel1, inverse=TRUE)
      if(!is.null(flowWorkspace::getTransformations(gating.set[[z]], channel=channel2, inverse=TRUE)))
      {
        c2.invTrans[[z]] = flowWorkspace::getTransformations(gating.set[[z]], channel=channel2, inverse=TRUE)
      }
      else  c2.invTrans[[z]] = flowWorkspace::getTransformations(gating.set[[z]], channel=channel2, inverse=TRUE)

      versions[z]<-as.numeric(flowWorkspace::getData(gating.set[[z]])@description$FCSversion)
      etim = as.POSIXlt(flowWorkspace::getData(gating.set[[z]])@description$`$ETIM`,format="%H:%M:%S")
      btim = as.POSIXlt(flowWorkspace::getData(gating.set[[z]])@description$`$BTIM`,format="%H:%M:%S")
      timediff[z] = as.numeric(difftime(etim, btim, units="min"))

      nodelist<-flowWorkspace::getNodes(gating.set[[z]], path = 1)
      node<-nodelist[which(nodelist==gate.name)]
      flowList[[z]] = flowWorkspace::getData(gating.set[[z]],node)

      if(use.beads)
      {

        bds<-nodelist[which(nodelist==beads.gate)]
        beadsList[[z]] = flowWorkspace::getData(gating.set[[z]],bds)

      }
    }

    flowWorkspace::closeWorkspace(workspace)


    for (p in 1:length(flowList)){

      if(invstate=="on")
      {

        if(use.beads)
        {


          while(i<2)
          {

            for (e in 1:length(beadsList)) {

              if(mean(versions)==2)
              {

                mean_chan1[e] = mean(log10(c1.invTrans[[e]](flowCore::exprs(beadsList[[e]])[,channel1])))
                mean_chan2[e] = mean(log10(c2.invTrans[[e]](flowCore::exprs(beadsList[[e]])[,channel2])))

              }

              else
              {

                mean_chan1[e] = mean(log10(c1.invTrans[[e]](flowCore::exprs(beadsList[[e]])[,channel1])/26.2144))
                mean_chan2[e] = mean(log10(c2.invTrans[[e]](flowCore::exprs(beadsList[[e]])[,channel2])/26.2144))

              }

            }

          }







          fmean_chan1[p] = mean_chan1[p]-mean(mean_chan1)
          fmean_chan2[p] = mean_chan2[p]-mean(mean_chan2)


          if(mean(versions)==2)
          {

            chan1 = log10(c1.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel1])) + fmean_chan1[p]
            chan2 = log10(c2.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel2])) + fmean_chan2[p]

          }

          else
          {

            chan1 = log10(c1.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel1])/26.144) + fmean_chan1[p]
            chan2 = log10(c2.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel2])/26.144) + fmean_chan2[p]

          }

        }
        else
        {
          if(mean(versions)==2)
          {
            chan1 = log10(c1.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel1]))
            chan2 = log10(c2.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel2]))
          }

          else

          {
            chan1 = log10(c1.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel1])/26.2144)
            chan2 = log10(c2.invTrans[[p]](flowCore::exprs(flowList[[p]])[,channel2])/26.2144)

          }

        }

        if(doplot){

          plot(chan2, chan1, pch=".", xlim=c(0,4), ylim=c(0,4), main=id[p], xlab=channel2, ylab=channel1, cex=3, xaxs="i",yaxs="i" )
          grid(nx=nbins, col="blue")

        }



        chan1.cuts <- seq(from = 0, to = 4, length = nbins + 1)
        chan2.cuts <- seq(from = 0, to = 4, length = nbins +1)
        index.chan1 <- cut(chan1, chan1.cuts, include.lowest = TRUE)
        index.chan2 <- cut(chan2, chan2.cuts, include.lowest = TRUE)
        m <- tapply(chan2, list(index.chan1, index.chan2), base::length)
        m[is.na(m)] <- 0
        matrixList[[p]] = m

      }

      if(invstate=="off")
      {

        if(use.beads)
        {


          while(i<2)
          {

            for (e in 1:length(beadsList)) {

              if(mean(versions)==2)
              {

                mean_chan1[e] = mean(log10(flowCore::exprs(beadsList[[e]])[,channel1]))
                mean_chan2[e] = mean(log10(flowCore::exprs(beadsList[[e]])[,channel2]))

              }

              else
              {

                mean_chan1[e] = mean(log10(flowCore::exprs(beadsList[[e]])[,channel1])/26.2144)
                mean_chan2[e] = mean(log10(flowCore::exprs(beadsList[[e]])[,channel2])/26.2144)

              }

            }

          }







          fmean_chan1[p] = mean_chan1[p]-mean(mean_chan1)
          fmean_chan2[p] = mean_chan2[p]-mean(mean_chan2)


          if(mean(versions)==2)
          {

            chan1 = log10(flowCore::exprs(flowList[[p]])[,channel1]) + fmean_chan1[p]
            chan2 = log10(flowCore::exprs(flowList[[p]])[,channel2]) + fmean_chan2[p]

          }

          else
          {

            chan1 = log10(flowCore::exprs(flowList[[p]])[,channel1]/26.144) + fmean_chan1[p]
            chan2 = log10(flowCore::exprs(flowList[[p]])[,channel2]/26.144) + fmean_chan2[p]

          }

        }
        else
        {
          if(mean(versions)==2)
          {
            chan1 = log10(flowCore::exprs(flowList[[p]])[,channel1])
            chan2 = log10(flowCore::exprs(flowList[[p]])[,channel2])
          }

          else

          {
            chan1 = log10(flowCore::exprs(flowList[[p]])[,channel1]/26.2144)
            chan2 = log10(flowCore::exprs(flowList[[p]])[,channel2]/26.2144)

          }

        }

        if(doplot){

          plot(chan2, chan1, pch=".", xlim=c(0,4), ylim=c(0,4), main=id[p], xlab=channel2, ylab=channel1, cex=3, xaxs="i",yaxs="i" )
          grid(nx=nbins, col="blue")

        }



        chan1.cuts <- seq(from = 0, to = 4, length = nbins + 1)
        chan2.cuts <- seq(from = 0, to = 4, length = nbins +1)
        index.chan1 <- cut(chan1, chan1.cuts, include.lowest = TRUE)
        index.chan2 <- cut(chan2, chan2.cuts, include.lowest = TRUE)
        m <- tapply(chan2, list(index.chan1, index.chan2), base::length)
        m[is.na(m)] <- 0
        matrixList[[p]] = m
      }


    }

    for (k in 1:length(matrixList))
    {
      matrixList[[k]] = as.matrix(as.data.frame(matrixList[[k]]))

    }

    MatList[[i]]=matrixList
    id.files = c(id.files,id)
    timediff.all = c(timediff.all,timediff)

  }



  for(y in 1:length(MatList))
  {

    for(w in 1:length(MatList[[y]]))

    {

      samples[[length(samples)+1]]= MatList[[y]][[w]]

    }

  }

  if(!is.null(flowrate)&&!is.null(volume)) volume=NULL

  if(is.null(flowrate)) vol1=rep(1, nsamples)
  else vol1=max(timediff.all*flowrate)/(timediff.all*flowrate)
  if (is.null(volume)) vol2=rep(1, nsamples)
  else vol2=max(volume)/volume



  for(d in 1:length(samples))

    samples[[d]] = vol1[d]*vol2[d]*dil[d]*samples[[d]]



  unmatrix <-function (x, byrow = FALSE)
  {
    rnames <- rownames(x)
    cnames <- colnames(x)
    if (is.null(rnames))
      rnames <- paste("r", 1:nrow(x), sep = "")
    if (is.null(cnames))
      cnames <- paste("c", 1:ncol(x), sep = "")
    nmat <- outer(rnames, cnames, paste, sep = ":")
    if (byrow) {
      vlist <- c(t(x))
      names(vlist) <- c(t(nmat))
    }
    else {
      vlist <- c(x)
      names(vlist) <- c(nmat)
    }
    return(vlist)
  }



  if (length(samples)==1)runbeta="FALSE"




  if (runalpha){
    alpha <- c()
    pielou <- c()
    for(w in 1:length(samples)){
      samps <- unmatrix(samples[[w]])
      alpha[w]<- vegan::diversity(samps, index=ialpha)
      pielou[w]<- vegan::diversity(samps)/log(vegan::specnumber(samps))
    }
  }


  if(runbeta){

    matriz <- matrix(ncol=length(samples), nrow=length(samples[[1]]))
    colnames(matriz) <- id.files

    for(g in 1:length(samples)) matriz[,g] <- unmatrix(samples[[g]], byrow="FALSE")


    beta = vegan::vegdist(t(matriz), method=ibeta)
    beta = as.matrix(beta)
  }


  indices <- list()

  if(runalpha){
    indices$alpha=alpha
    names(indices$alpha)=id.files
    indices$pielou=pielou
    names(indices$pielou)=id.files

  }


  if(runbeta){
    indices$beta=beta

  }

  if(save.csv) write.csv(indices, file=file.name, row.names=TRUE)

  return(indices)

}
