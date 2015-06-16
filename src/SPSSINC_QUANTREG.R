#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 1989, 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# version 1.2.3

# history
# 28-sep-2013 add final newline to file
# 07-jul-2014 html help

helptext="The SPSSINC QUANTREG command requires the R Integration Plug-in 
and the R quantreg package.

SPSSINC QUANTREG DEPENDENT=dependent variable 
ENTER=independent variables [QUANTILES=list of quantiles]
[/OPTIONS [MISSING={LISTWISE**}]
                   {FAIL      }
[METHOD={BR**}]
        {FN  }
        {PFN }
[STDERR={RANK**}]
        {IID   }
        {NID   }
        {KER   }
        {BOOT  }
[PLOT]
[EXECUTE={TRUE**}] ]
         {FALSE }
[/SAVE [RESIDUALSDATASET=datasetname] [COEFSDATASET=datasetname] 
       [PROGRAMFILE=filespec] ]

Split files and weight are not honored by this command.

SPSSINC QUANTREG /HELP prints this information and does nothing else.

Example:
SPSSINC QUANTREG DEPENDENT=mpg ENTER=engine weight.

Execute the rq function for quantile regression from the R quantreg package.
DEPENDENT and ENTER specify the dependent and independent
variable names.  

Categorical independent variables are automatically converted
appropriately to factors.  A constant term is automatically included.

QUANTILES is a list of quantiles with values between 0 and 1.  The quantile
regression is performed for each value.  The default is .5

MISSING=LISTWISE causes listwise deletion of missing values (but
missings in factors are included as a level).  FAIL stops the procedure if
missing values are encountered.

METHOD specifies the computational algorithm.
BR is Barrodale and Roberts, recommended for up to several thousand cases.
FN is Frisch-Newton, recommended for large problems.
PFN is Frisch-Newton after preprocessing, recommended for very large problems.

STDERR specifies how to compute coefficient standard errors.
RANK specifies confidence intervals by inverting a rank test
IID assumes iid errors and computes asymptotic estimates
NID computes a Huber sandwich estimate
KER uses a kernel estimate
BOOT computes a bootstrapped estimator

Details can be found in the R help for rq.

PLOT causes each coefficient value to be plotted against the specified quantile
values if there is more than one quantile.

EXECUTE=FALSE runs the command syntax without running the quantile regression.  
This is mainly useful in combination with SAVE PROGRAMFILE.

/SAVE RESIDUALSDATASET causes a  dataset containing the residuals to be created.
The name must not already be in use.
The case number is included as cases will only be written for input cases with no
missing data.  A variable is created for each regression.

COEFSDATASET causes a dataset containing the coefficients to be created.
The name must not already be in use.
A variable is created for each regression.

PROGRAMFILE causes the R code that implements the quantile regression to be 
written to the specified file.  Since the rq function has features not exposed
in this extension command, the generated program can be a useful starting point 
for additional specifications.
"

quantreg<-function(dep, enter, taus=.5, missing="listwise", method="br", stderr="rank",
            plot=FALSE, residualsdataset=NULL, coefsdataset=NULL){

    domain<-"SPSSINC_QUANTREG"
    setuplocalization(domain)
    
    tryCatch(library(quantreg), error=function(e){
        stop(gettextf("The R %s package is required but could not be loaded.","quantreg",domain=domain),call.=FALSE)
        }
    )

    if (identical(missing,"listwise")) {missing<-na.exclude} else {missing<-na.fail}
    
    allvars <- c(dep,enter)
    model <- paste(dep,"~",paste(enter,collapse="+"))

    dta<-spssdata.GetDataFromSPSS(allvars,missingValueToNA=TRUE,factorMode="labels")
    
    res <- tryCatch(
            summary(rqres<-rq(as.formula(model),data=dta, tau=taus, 
              na.action=missing, method=method), se=stderr),
            error=function(e) {return(c(gettext("ERROR:",domain=domain),e))}
    )

    if (!is.null(res$message)) {print(res$message)} else {
        miss<-ifelse(identical(missing,na.exclude),"na.exclude","na.fail")
        StartProcedure(gettext("Quantile Regression",domain=domain),"SPSSINC QUANTREG")
        for (i in 1:length(taus)) {
            if (length(taus)>1) coeff<-res[[i]]$coefficients else coeff<-res$coefficients
            for (j in 1:length(attributes(coeff)$dimnames[[1]])){
                attributes(coeff)$dimnames[[1]][[j]]=gettext(attributes(coeff)$dimnames[[1]][[j]],domain=domain)
            }
            for (j in 1:length(attributes(coeff)$dimnames[[2]])){
                attributes(coeff)$dimnames[[2]][[j]]=gettext(attributes(coeff)$dimnames[[2]][[j]],domain=domain)
            }
            spsspivottable.Display(coeff, 
                title=gettextf("Coefficients, quantile = %s",taus[i],domain=domain),templateName="SpssincQuantreg",
                caption=paste("rq(formula = ",model,", tau = ",deparse(taus),", data = dta, na.action = ",miss,
                    ", method = ",dQuote(method),")",sep=""),isSplit=FALSE)
        }
        spsspkg.EndProcedure()
    }

    if (is.null(rqres$message)) {
        if (!is.null(residualsdataset)){
            residualsdict = spssdictionary.CreateSPSSDictionary(c("caseNumber", gettext("Case Number",domain=domain), 0, "F8.0", "nominal"))
            for (i in 1:length(taus)){
                residualsdict<-spssdictionary.CreateSPSSDictionary(residualsdict,
                    c(paste("Residuals_",i-1,sep=""), gettextf("tau = %s",taus[i],domain=domain), 0 ,"F8.2", "scale"))
            }
            tryCatch({
                spssdictionary.SetDictionaryToSPSS(residualsdataset, residualsdict)
                resids = list()
                if (length(taus)>1)
                    for (i in 1:length(taus)) {resids[[i]] = residuals(rqres)[,i]}
                else 
                    resids[[1]] = residuals(rqres)
                resids = data.frame(resids)
                resids[is.na(resids)] <- NaN
                spssdata.SetDataToSPSS(residualsdataset, data.frame(row.names(resids), resids))
                }, 
                error=function(e) {print(e)
                cat(gettext("Failed to create residuals dataset. Dataset name must not already exist: ",domain=domain),residualsdataset)
                }
            )
        }
        if (!is.null(coefsdataset)){
            coefsdict = spssdictionary.CreateSPSSDictionary(c("term", gettext("Variable or Factor Level",domain=domain), 100, "A100", "nominal"))
            for (i in 1:length(taus)){
                coefsdict<-spssdictionary.CreateSPSSDictionary(coefsdict,
                    c(paste("Coefficients_",i-1,sep=""), gettextf("tau = %s",taus[i],domain=domain), 0, "F8.2", "scale"))
            }
            tryCatch({
                spssdictionary.SetDictionaryToSPSS(coefsdataset, coefsdict)
                coefs= list()
                if (length(taus)>1)
                    for (i in 1:length(taus)) {coefs[[i]] = coefficients(res[[i]])[,1]}
                else
                    coefs[[1]] = coefficients(res)[,1]
                coefs = data.frame(coefs)
                spssdata.SetDataToSPSS(coefsdataset, data.frame(row.names(coefs), coefs))
                }, 
                error=function(e) {print(e)
                cat(gettext("Failed to create coefficients dataset. Dataset name must not already exist: ",domain=domain),residualsdataset)
                }
            )

        }
            
        if (plot & length(taus) > 1) plot(res,ask=FALSE)
            
        spssdictionary.EndDataStep()
    }
    
    res <- tryCatch(rm(list=ls()),warning=function(e){return(NULL)})

}

StartProcedure<-function(procname, omsid){
if (as.integer(substr(spsspkg.GetSPSSVersion(),1, 2)) >= 19)
   spsspkg.StartProcedure(procname,omsid)
else
   spsspkg.StartProcedure(omsid)
}

caller<-function(dep, enter, taus=.5, missing="listwise",  method="br", stderr="rank", 
            plot=FALSE, residualsdataset=NULL, coefsdataset=NULL, programfile=NULL, execute=TRUE){

    if (length(taus) == 1)
        title<-"# SPSSINC QUANTREG, one quantile value\n"
    else
        title<-"# SPSSINC QUANTREG, multiple quantile values\n"
   
    if(!is.null(programfile)){
        lines<-c(title,
            "quantreg<-",
            attr(quantreg,"source"),
            paste("dep<-",dQuote(dep),sep=""),
            paste("enter<-",deparse(enter),sep=""),
            paste("taus<-",deparse(taus),sep=""),
            paste("missing<-",dQuote(missing),sep=""),
            paste("method<-",dQuote(method),sep=""),
            paste("stderr<-",dQuote(stderr),sep=""),
            paste("plot<-",plot,sep=""))
        func<-"quantreg(dep, enter, taus, missing, method, stderr, plot"
        if(!is.null(residualsdataset)){
            func<-paste(func,", residualsdataset=",dQuote(residualsdataset),sep="")
        }
        if(!is.null(coefsdataset)){
            func<-paste(func,", coefsdataset=",dQuote(coefsdataset),sep="")
        }
        func<-paste(func,")",sep="")
        lines<-c(lines,func)        
        f<-file(description=programfile,open="wb",encoding="UTF-8")
        writeLines(lines,con=f)
        close(f)
    }
    
    if (execute) quantreg(dep, enter, taus, missing, method, stderr, plot, residualsdataset, coefsdataset)
    
}

setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

Run<-function(args){
    cmdname = args[[1]]
    args <- args[[2]]
    oobj = spsspkg.Syntax(templ=list(
                spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep", islist=FALSE),
                spsspkg.Template("ENTER", subc="",  ktype="existingvarlist", var="enter", islist=TRUE),
                spsspkg.Template("QUANTILES", subc="", ktype="float", var="taus", islist=TRUE),
                spsspkg.Template("MISSING", subc="OPTIONS",ktype="str", var="missing"),
                spsspkg.Template("METHOD", subc="OPTIONS", ktype="str", var="method", islist=FALSE),
                spsspkg.Template("STDERR", subc="OPTIONS", ktype="str", var="stderr", islist=FALSE),
                spsspkg.Template("PLOT", subc="OPTIONS", ktype="bool", var="plot", islist=FALSE),
                spsspkg.Template("EXECUTE", subc="OPTIONS", ktype="bool", var="execute"),
                spsspkg.Template("PROGRAMFILE", subc="SAVE", ktype="literal", var="programfile"),
                spsspkg.Template("RESIDUALSDATASET", subc="SAVE", ktype="literal", var="residualsdataset"),
                spsspkg.Template("COEFSDATASET", subc="SAVE", ktype="literal", var="coefsdataset")
                ))

    if ("HELP" %in% attr(args,"names"))
        #writeLines(helptext)
        helper(cmdname)
    else
        res <- spsspkg.processcmd(oobj,args,"caller")
}


helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}