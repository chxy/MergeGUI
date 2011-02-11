##' obtain the intersection of a list
##'
##'
##' @param vname a
##' @param simplifiedname a
##' @return The outputs are 'public', 'individual', 'uniq', and
##' 'simpleuniq'.  'public' is a vector of the intersection of
##' 'simplifiedname'.  'individual' is a matrix with the original
##' colnames matched to 'public' in all files.  'simpleuniq' is a list
##' of the left part of 'simplifiedname' if we pick 'public' out.
##' 'uniq' is a list of the left part of 'vname' if we pick
##' 'individual' out.
##' @author Xiaoyue Cheng
##' @export
##' @examples ## TODO
intersect2 = function(vname, simplifiedname) {
    #####---------------------------------------------------------------#####
    ##  intersect2 is a function to .     ##
    ##  The inputs are both vname and simplifiedname, which are matched.   ##
    #####---------------------------------------------------------------#####

        s = as.vector(simplifiedname[[1]])
        for (i in 2:length(simplifiedname)) {
            s = intersect(as.vector(simplifiedname[[i]]), s)
        }
        v1 = matrix(nrow = length(s), ncol = length(vname))
        v2 = vname
        v3 = simplifiedname
        if (length(s) > 0) {
            for (i in 1:length(vname)) {
                for (j in 1:length(s)) {
                  if (s[j] %in% simplifiedname[[i]]) {
                    v1[j, i] = vname[[i]][which(simplifiedname[[i]] ==
                      s[j])[1]]
                    tmp = vname[[i]][-which(simplifiedname[[i]] ==
                      s[j])[1]]
                    v2[[i]] = intersect(v2[[i]], tmp)
                    v3[[i]] = v3[[i]][-which(v3[[i]] == s[j])[1]]
                  }
                }
            }
        }
        return(list(public = s, individual = v1, uniq = v2, simpleuniq = v3))
    }

##' Function to short the names from a template.
##'
##'
##' @param namevector a
##' @return varclass
##' @author Xiaoyue Cheng
##' @export
##' @examples ## TODO
simplifynames=function(namevector) {
		n=max(nchar(namevector))
		for (i in 1:n){
			if (!all(substr(namevector,i,i)==substr(namevector,i,i)[1])){
				newnamevec=substring(namevector,i)
				return(newnamevec)
			}
		}
		return(namevector)
	}

##' a function to detect the type of the variables.
##'
##'
##' @param nametable.class a
##' @param dataset.class a
##' @return a vector matching the rows of 'nametable'. The value
##' includes NA if any variable are only NA's.
##' @author Xiaoyue Cheng
##' @export
##' @examples ## TODO
var.class = function(nametable.class, dataset.class) {
    #####----------------------------------------------------------#####
    ##  var.class is   ##
	##  The value of this function is           ##
    #####----------------------------------------------------------#####
        varclass = rep("NA", nrow(nametable.class))
        for (i in 1:nrow(nametable.class)) {
            notNAcolset = which(!is.na(nametable.class[i, ]))
            notNAcol = NA
            for (k in notNAcolset) {
                if (sum(!is.na(dataset.class[[k]][, nametable.class[i,
                  k]])) > 0)
                  notNAcol = k
            }
            if (!is.na(notNAcol)) {
                varclass[i] = class(dataset.class[[notNAcol]][,
                  nametable.class[i, notNAcol]])
            }
        }
        return(varclass)
    }


##' function
##'
##'
##' @param nametable.class a
##' @param dataset.class a
##' @param name.class a
##' @return NULL
##' @author Xiaoyue Cheng
##' @export
##' @examples ## TODO
scale.rpart = function(nametable.class, dataset.class, name.class) {
		txtpb = txtProgressBar(min=0,max=1,width = 100,style=3)
	    varclass = var.class(nametable.class,dataset.class)
		rows = unlist(lapply(dataset.class, nrow))
		selectedvariables = which(varclass %in% c('numeric','integer','logical','factor'))
		mergedata = matrix(nrow = sum(rows), ncol = length(selectedvariables) + 1)
        colnames(mergedata) = c("source", name.class[selectedvariables])
		setTxtProgressBar(txtpb, 0.1)
        for (i in 1:length(dataset.class)) {
            tmp = matrix(c(rep(colnames(nametable.class)[i], rows[i]),
                rep(NA, rows[i] * length(selectedvariables))), nrow = rows[i])
            colnames(tmp) = c("source", nametable.class[selectedvariables, i])
            tmp[, na.omit(nametable.class[selectedvariables, i])] = as.matrix(dataset.class[[i]])[,
                na.omit(nametable.class[selectedvariables, i])]
            mergedata[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i],
                ] = tmp
        }
		setTxtProgressBar(txtpb, 0.2)
		mergedata=as.data.frame(mergedata)
		mergedata$source=factor(mergedata$source)
		num = which(varclass[selectedvariables] %in% c('integer','numeric'))
		fac = which(varclass[selectedvariables] %in% c('logical','factor'))
		if (length(num)!=0) {
			mergedata[,num+1] = sapply(mergedata[,num+1],function(avec){
				as.numeric(as.character(avec))})
		}
		if (length(fac)!=0) {
			mergedata[,fac+1]=sapply(mergedata[,fac+1],as.factor)
		}
		setTxtProgressBar(txtpb, 0.3)
		#fit_rpart = rpart(source~., data=mergedata, method="class",
		#	control=rpart.control(maxcompete=ncol(mergedata)-1))
		#rpartout = rownames(as.data.frame(fit_rpart$splits))
		#scaleclass = round(as.data.frame(fit_rpart$splits)[name.class,'improve'],2)
		#scaleclass[is.na(scaleclass)] = ''
		#return(scaleclass)

		res = rep(9, nrow(nametable.class))
		for (i in 2:ncol(mergedata)) {
			fit_rpart = rpart(mergedata$source~mergedata[,i], control=c(maxdepth=1))
                        group = mergedata$source
			tmperror = weighted.mean(residuals(fit_rpart), 1/table(group)[group])
			res[name.class==colnames(mergedata)[i]] = round(tmperror,3)
			if (tmperror==0){
				if (all(is.na(mergedata[fit_rpart$where==2,i])) |
					all(is.na(mergedata[fit_rpart$where==3,i]))) {
					res[name.class==colnames(mergedata)[i]] = 9
				}
			}
			setTxtProgressBar(txtpb, 0.3+0.65*i/ncol(mergedata))
		}
		setTxtProgressBar(txtpb, 1)
		return(as.character(res))
            }

##' function
##'
##'
##' @param nametable.class a
##' @param dataset.class a
##' @param name.class a
##' @return NULL
##' @author Xiaoyue Cheng
##' @export
##' @examples ## TODO
scale.kstest = function(nametable.class, dataset.class, name.class) {
		txtpb = txtProgressBar(min=0,max=1,width = 100,style=3)
	    varclass = var.class(nametable.class,dataset.class)
		rows = unlist(lapply(dataset.class, nrow))
		selectedvariables = which(varclass %in% c('numeric','integer'))
		mergedata = matrix(nrow = sum(rows), ncol = length(selectedvariables) + 1)
        colnames(mergedata) = c("source", name.class[selectedvariables])
		setTxtProgressBar(txtpb, 0.05)
        for (i in 1:length(dataset.class)) {
            tmp = matrix(c(rep(colnames(nametable.class)[i], rows[i]),
                rep(NA, rows[i] * length(selectedvariables))), nrow = rows[i])
            colnames(tmp) = c("source", nametable.class[selectedvariables, i])
            tmp[, na.omit(nametable.class[selectedvariables, i])] = as.matrix(dataset.class[[i]])[,
                na.omit(nametable.class[selectedvariables, i])]
            mergedata[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i],
                ] = tmp
        }
		setTxtProgressBar(txtpb, 0.1)
		mergedata=as.data.frame(mergedata)
		mergedata$source=factor(mergedata$source)
        mergedata[,2:ncol(mergedata)]=sapply(mergedata[,2:ncol(mergedata)],
			function(avec){as.numeric(as.character(avec))})
		setTxtProgressBar(txtpb, 0.15)

		kstestout = c()
		scaleclass = rep(9, nrow(nametable.class))
		for (i in 2:ncol(mergedata)) {
			tmpdat = mergedata[,c(1,i)]
			sig = c()
			for (j in 1:length(levels(tmpdat[,1]))) {
				a = scale(tmpdat[tmpdat[,1]==levels(tmpdat[,1])[j],2], scale=FALSE)
				b = scale(tmpdat[tmpdat[,1]!=levels(tmpdat[,1])[j],2], scale=FALSE)
				if (all(is.na(a)) | all(is.na(b))) {
					sig[j] = 1
				} else {
					sig[j] = ks.test(a, b)$p.value
				}
			}
			if (!all(sig>0.1)) {
				kstestout = c(kstestout,i)
				scaleclass[selectedvariables][i-1] = min(sig, na.rm=TRUE)
			}
			setTxtProgressBar(txtpb, 0.15+0.8*i/ncol(mergedata))
		}
		setTxtProgressBar(txtpb, 1)
		return(as.character(round(scaleclass,3)))
	}
##' function
##'
##'
##' @param nametable.class a
##' @param dataset.class a
##' @param name.class a
##' @return NULL
##' @author Xiaoyue Cheng
##' @export
##' @examples ## TODO
scale.missing = function(nametable.class, dataset.class, name.class) {
		txtpb = txtProgressBar(min=0,max=1,width = 100,style=3)
		rows = unlist(lapply(dataset.class, nrow))
		mergedata = matrix(nrow = sum(rows), ncol = length(name.class) + 1)
        colnames(mergedata) = c("source", name.class)
        for (i in 1:length(dataset.class)) {
            tmp = matrix(c(rep(colnames(nametable.class)[i], rows[i]),
                rep(NA, rows[i] * length(name.class))), nrow = rows[i])
            colnames(tmp) = c("source", nametable.class[, i])
            tmp[, na.omit(nametable.class[, i])] = as.matrix(dataset.class[[i]])[,
                na.omit(nametable.class[, i])]
            mergedata[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i],
                ] = tmp
        }
		mergedata=as.data.frame(mergedata)
		mergedata$source=factor(mergedata$source)
		setTxtProgressBar(txtpb, 0.1)

		chi2testout = c()
		missingclass = rep(9, nrow(nametable.class))
		for (i in 2:ncol(mergedata)) {
			tmpdat = mergedata[,c(1,i)]
			missingcount = matrix(0, nrow=length(levels(tmpdat[,1])), ncol=2)
			for (j in 1:length(levels(tmpdat[,1]))) {
				missingcount[j,1] = sum(is.na(tmpdat[tmpdat[,1]==levels(tmpdat[,1])[j],2]))
				missingcount[j,2] = rows[j] - missingcount[j,1]
			}
			missingclass[i-1] = chisq.test(missingcount)$p.value
			setTxtProgressBar(txtpb, 0.1+0.8*i/ncol(mergedata))
		}
		setTxtProgressBar(txtpb, 1)
		return(as.character(round(missingclass,3)))
	}

##' function
##'
##'
##' @return NULL
##' @author Xiaoyue Cheng
##' @export
##' @examples ## TODO
mergeGUI = function() {
	mergegui_env = new.env()

mergefunc = function(h, ...) {

    undo = function(h, ...) {
    #####-----------------------------------------------------#####
	##  The following buttons are used for switching variables.  ##
    ##  undo button.                                             ##
    #####-----------------------------------------------------#####
        mergegui_env$idx <- mergegui_env$idx - 1
        if (mergegui_env$idx == 0) {
            gmessage("You can not undo anymore!")
            mergegui_env$idx <- 1
        }
        for (i in 1:n) {
            gt2[[i]][] = mergegui_env$hstry1[[mergegui_env$idx]][, i]
        }
        mergegui_env$redo.indicate <- 1
        gt4[, 2] = mergegui_env$hstry2[[mergegui_env$idx]]
        gt4[, 3] = mergegui_env$hstry3[[mergegui_env$idx]]
        gt5[, 2] <- gt4[, 2]
        gt5[, 3] <- gt4[, 3]
    }

    redo = function(h, ...) {
    #####-----------------------------------------------------#####
	##  The following buttons are used for switching variables.  ##
    ##  redo button.                                             ##
    #####-----------------------------------------------------#####
        if (mergegui_env$redo.indicate == 0) {
            gmessage("There is nothing to redo.")
            return()
        }
        mergegui_env$idx <- mergegui_env$idx + 1
        if (mergegui_env$idx > length(mergegui_env$hstry1)) {
            gmessage("You can not redo anymore!")
            mergegui_env$idx <- length(mergegui_env$hstry1)
        }
        for (i in 1:n) {
            gt2[[i]][] = mergegui_env$hstry1[[mergegui_env$idx]][, i]
        }
        gt4[, 2] = mergegui_env$hstry2[[mergegui_env$idx]]
        gt4[, 3] = mergegui_env$hstry3[[mergegui_env$idx]]
        gt5[, 2] <- gt4[, 2]
        gt5[, 3] <- gt4[, 3]
    }

    reset = function(h, ...) {
    #####-----------------------------------------------------#####
	##  The following buttons are used for switching variables.  ##
    ##  reset button.                                            ##
    #####-----------------------------------------------------#####
        for (i in 1:n) {
            gt2[[i]][] = mergegui_env$hstry1[[1]][, i]
        }
        mergegui_env$redo.indicate = 1
        gt4[, 2] = nameintersect
        gt4[, 3] = var.class(nametable,dataset)
        gt5[, 2] <- gt4[, 2]
        gt5[, 3] <- gt4[, 3]
    }

	VariableOptions = function(h, ...) {
    #####------------------------------------------------------#####
    ##  VariableOptions is the handler when double clicking gt4.  ##
	##  It gives a new window for                                 ##
	##          editing the attributes of variables.              ##
    #####------------------------------------------------------#####
        gt4input0 = gwindow("Attributes", visible = T, width = 300,
            height = 200)
        gt4input = ggroup(horizontal = FALSE, container = gt4input0,
            expand = TRUE)
        gt4input1 = ggroup(container = gt4input, expand = TRUE)
        gt4input2 = ggroup(container = gt4input, expand = TRUE)
		gt4input4 = ggroup(container = gt4input, horizontal = FALSE, expand = TRUE)
        gt4input3 = ggroup(container = gt4input, expand = TRUE)

        gt4input11 = glabel("Name:", container = gt4input1)
        gt4input12 = gedit(text = svalue(gt4), container = gt4input1,
            expand = TRUE)
        gt4input21 = glabel("Class:", container = gt4input2)
        gt4input22 = gcombobox(union(gt4[svalue(gt4, index = TRUE),
            3], c("integer", "numeric", "character", "factor")),
            container = gt4input2, expand = TRUE, handler = function(h,...){
			# if (svalue(gt4input22) %in% c('character','factor')) {
				# svalue(gt4input41) = FALSE
				# enabled(gt4input41) = FALSE
				# enabled(gt4input42) = FALSE
				# for (i in 1:n) {
					# svalue(gt4input42edit[[i]])=""
				# }
			# } else {
				# enabled(gt4input41) = TRUE
			#}
			})

		# gt4input41 = gcheckbox("Different scales", checked =
			# (gt4[svalue(gt4, index = TRUE),4]=='X'), container = gt4input4,
			# expand = TRUE, handler = function(h,...){
			# if (svalue(gt4input22) %in% c('integer','numeric')) {
				# enabled(gt4input42) = svalue(gt4input41)
				# if (!svalue(gt4input41)) {
					# for (i in 1:n) {
						# svalue(gt4input42edit[[i]])=""
					# }
				# }
			# }
			# })
		# gt4input42 = gframe("Formulas to unify the scales",
			# horizontal = FALSE, container =gt4input4,
			# expand=TRUE)
		# enabled(gt4input42) = svalue(gt4input41)
		# gt4input42group = gt4input42label = gt4input42edit = list()
		# for (i in 1:n) {
			# gt4input42group[[i]]=ggroup(container = gt4input42, expand = TRUE)
			# gt4input42label[[i]]=glabel(sub('.csv','', basename(gtfile)[i]),
				# container = gt4input42group[[i]])
			# gt4input42edit[[i]]=gedit(expand = TRUE, container = gt4input42group[[i]])
		# }

        gt4input31 = gbutton("Ok", container = gt4input3, expand = TRUE,
            handler = function(h, ...) {
                if (svalue(gt4input12) != "") {
                  gt4[svalue(gt4, index = TRUE), 2] = svalue(gt4input12)
                  gt4[svalue(gt4, index = TRUE), 3] = svalue(gt4input22)
				  # if (svalue(gt4input41)) {
					# gt4[svalue(gt4, index = TRUE), 4] = "X"
				  # } else {
				    # gt4[svalue(gt4, index = TRUE), 4] = ""
				  # }
                  gt5[, 2] = gt4[, 2]
                  gt5[, 3] = gt4[, 3]
                  mergegui_env$hstry2[[mergegui_env$idx]] <- gt4[, 2]
                  mergegui_env$hstry3[[mergegui_env$idx]] <- gt4[, 3]
                  dispose(gt4input0)
                }
                else {
                  gmessage("Variable name could not be empty!")
                }
            })
        gt4input32 = gbutton("Cancel", container = gt4input3,
            expand = TRUE, handler = function(h, ...) {
                dispose(gt4input0)
            })
    }

    smmry = function(h, ...) {
	#####---------------------------------#####
    ##  smmry is the handler of gbcombo431.  ##
	##  (gbutton: Numeric Summary)           ##
    #####---------------------------------#####
        graphics.off()
        name.select = svalue(gt4, index = TRUE)

        if (length(name.select) == 0) {
            gmessage("Please select the variables!")
            return()
        }
        name.table = matrix(nrow = length(name.select), ncol = n)
        for (i in 1:n) {
            name.table[, i] = gt2[[i]][gt4[name.select,1], ]
        }
        name.intersect = as.vector(svalue(gt4))
        name.class = gt4[name.select, 3]
        summarytable = list()
        for (i in 1:length(name.select)) {
            if (name.class[i] != "NA") {
                if (name.class[i] == "numeric" | name.class[i] ==
                  "integer") {
                  summarytable[[i]] = matrix(NA, ncol = n, nrow = 7,
                    dimnames = list(c("size", "NA#s", "mean",
                      "std", "min", "median", "max"),
					  simplifynames(gsub('.csv','',basename(gtfile)))))
                  names(summarytable)[i] = name.intersect[i]
                  for (j in 1:n) {
                    if (!is.na(name.table[i, j])) {
                      tmpdata = dataset[[j]][, name.table[i,
                        j]]
                      summarytable[[i]][1, j] = length(tmpdata)
                      summarytable[[i]][2, j] = sum(is.na(tmpdata))
                      summarytable[[i]][3, j] = mean(tmpdata,
                        na.rm = TRUE)
                      summarytable[[i]][4, j] = sd(tmpdata, na.rm = TRUE)
                      summarytable[[i]][5, j] = min(tmpdata,
                        na.rm = TRUE)
                      summarytable[[i]][6, j] = median(tmpdata,
                        na.rm = TRUE)
                      summarytable[[i]][7, j] = max(tmpdata,
                        na.rm = TRUE)
                    }
                  }
				  summarytable[[i]] = data.frame(t(summarytable[[i]]))
				  summarytable[[i]][,1] = as.integer(as.character(summarytable[[i]][,1]))
				  summarytable[[i]][,2] = as.integer(as.character(summarytable[[i]][,2]))
				  summarytable[[i]][,3] = as.character(round(summarytable[[i]][,3]),3)
				  summarytable[[i]][,4] = as.character(round(summarytable[[i]][,4]),3)
				  summarytable[[i]][,5] = as.character(round(summarytable[[i]][,5]),3)
				  summarytable[[i]][,6] = as.character(round(summarytable[[i]][,6]),3)
				  summarytable[[i]][,7] = as.character(round(summarytable[[i]][,7]),3)
				  summarytable[[i]] = cbind(File=rownames(summarytable[[i]]),summarytable[[i]])
                }
                else {
                  summarytable[[i]] = matrix(NA, ncol = n, nrow = 10,
                    dimnames = list(c("size", "NA#s", "levels",
                      "matched levels", "top 1 level", "amount 1",
                      "top 2 level", "amount 2", "top 3 level", "amount 3"),
                      simplifynames(gsub('.csv','',basename(gtfile)))))
                  names(summarytable)[i] = name.intersect[i]
                  matchedlevels = list()
                  for (j in 1:n) {
                    if (!is.na(name.table[i, j])) {
                      if (sum(!is.na(dataset[[j]][, name.table[i,
                        j]])) > 0) {
                        matchedlevels[[j]] = names(table(dataset[[j]][,
                          name.table[i, j]], useNA = "no"))
                      }
                      else {
                        matchedlevels[[j]] = NA
                      }
                    }
                    else {
                      matchedlevels[[j]] = NA
                    }
                  }
                  mtch = intersect2(matchedlevels, matchedlevels)
                  for (j in 1:n) {
                    if (!is.na(name.table[i, j])) {
                      tmpdata = dataset[[j]][, name.table[i,
                        j]]
                      tmptable = sort(table(tmpdata, useNA = "no"),
                        decreasing = TRUE)
                      summarytable[[i]][1, j] = length(tmpdata)
                      summarytable[[i]][2, j] = sum(is.na(tmpdata))
                      summarytable[[i]][3, j] = length(tmptable)
                      summarytable[[i]][4, j] = length(mtch$public)
                      summarytable[[i]][5, j] = names(tmptable)[1]
                      summarytable[[i]][6, j] = tmptable[1]
                      if (length(tmptable) > 1) {
                        summarytable[[i]][7, j] = names(tmptable)[2]
                        summarytable[[i]][8, j] = tmptable[2]
                      }
                      if (length(tmptable) > 2) {
                        summarytable[[i]][9, j] = names(tmptable)[3]
                        summarytable[[i]][10, j] = tmptable[3]
                      }
                    }
                  }
				  summarytable[[i]] = data.frame(t(summarytable[[i]]))
				  summarytable[[i]][,1] = as.integer(as.character(summarytable[[i]][,1]))
				  summarytable[[i]][,2] = as.integer(as.character(summarytable[[i]][,2]))
				  summarytable[[i]][,3] = as.integer(as.character(summarytable[[i]][,3]))
				  summarytable[[i]][,4] = as.integer(as.character(summarytable[[i]][,4]))
				  summarytable[[i]][,6] = as.integer(as.character(summarytable[[i]][,6]))
				  summarytable[[i]][,8] = as.integer(as.character(summarytable[[i]][,8]))
				  summarytable[[i]][,10] = as.integer(as.character(summarytable[[i]][,10]))
				  summarytable[[i]] = cbind(File=rownames(summarytable[[i]]),summarytable[[i]])
                }
            }
            else {
                summarytable[[i]] = matrix(NA, nrow = n, ncol = 1,
                  dimnames = list(simplifynames(gsub('.csv','',basename(gtfile))), NULL))
                names(summarytable)[i] = name.intersect[i]
            }
        }

		gbcombo441 = list()

			delete(group43, group45)
			group45 <- ggroup(cont=group43,expand = TRUE, use.scrollwindow = TRUE)
			gbcombo44 <- glayout(container = group45,expand = TRUE, use.scrollwindow = TRUE)

		for (i in 1:length(name.select)){
			gbcombo44[i*2-1, 1] = glabel(names(summarytable)[i],cont=gbcombo44)
			gbcombo44[i*2, 1, expand = TRUE] = gbcombo441[[i]] = gtable(summarytable[[i]], container = gbcombo44)
		}

    }

    graph = function(h, ...) {
    #####---------------------------------#####
    ##  graph is the handler of gbcombo432.  ##
	##  (gbutton: Graphic Summary)           ##
    #####---------------------------------#####
        graphics.off()
			delete(group43, group45)
			group45 <- ggroup(cont=group43,expand = TRUE, use.scrollwindow = TRUE)
			gbcombo44 <- glayout(container = group45,expand = TRUE, use.scrollwindow = TRUE)
		gbcombo44[1, 1, expand = TRUE] = gbcombo442 = ggroup(container = gbcombo44, use.scrollwindow = TRUE)

		yscale = svalue(radio121)
        name.select = svalue(gt4, index = TRUE)
        if (length(name.select) == 0) {
            gmessage("Please select one variables!")
            return()
        }
        if (length(name.select) > 1) {
            gmessage("You selected more than one variable!
				Only the first two numeric variables are shown.")
			whichnumeric = gt4[name.select,3] %in% c('integer','numeric')
			if (!any(whichnumeric)) {
				name.select = name.select[1]
			} else {
				if (sum(whichnumeric) > 1) {
					name.select = name.select[whichnumeric][1:2]
				} else {
					name.select = name.select[whichnumeric][1]
				}
			}
        }

		if (length(name.select)==1) {
        name.table = rep(NA, n)
        for (i in 1:n) {
            name.table[i] = gt2[[i]][gt4[name.select,1], ]
        }
        name.intersect = as.character(gt4[name.select,2])
        name.class = as.character(gt4[name.select, 3])
        mergedata = data.frame(source = rep(simplifynames(gsub('.csv','',basename(gtfile))),
            rows))

        is.num = FALSE
        if (name.class != "NA") {
            if (name.class == "numeric" | name.class == "integer") {
                is.num = TRUE
                tmp.num = c()
                for (i in 1:n) {
                  if (!is.na(name.table[i])) {
                    tmp.num = c(tmp.num, dataset[[i]][, name.table[i]])
                  }
                  else {
                    tmp.num = c(tmp.num, rep(NA, rows[i]))
                  }
                }
                mergedata = data.frame(mergedata, as.numeric(tmp.num))
				mergedata[,1] = reorder(mergedata[,1], mergedata[,2], median, na.rm=TRUE)
            }
            else {
                tmp.chr = rep("na", sum(rows))
                for (i in 1:n) {
                  if (!is.na(name.table[i])) {
                    tmp.chr[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i]] = as.character(dataset[[i]][,
                      name.table[i]])
                  }
                  else {
                    tmp.chr[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i]] = rep(NA,
                      rows[i])
                  }
                }
                levelorder = names(sort(table(tmp.chr), decreasing = FALSE))
                mergedata = cbind(mergedata, factor(tmp.chr,
                  levels = levelorder))
            }
        }
        else {
            mergedata = cbind(mergedata, rep(NA, sum(rows)))
        }


		gbcombo4421 = ggraphics(container = gbcombo442,
            height = ifelse(is.num, 75 * 3 * n, 75 * 6),  expand = TRUE)

        colnames(mergedata) = c("source", name.intersect)

		if (yscale=="regular y scale") {
			eval(parse(text = paste("print(qplot(", name.intersect,
				",data=mergedata,facets=", ifelse(is.num,
				"source~.)", "~source)+coord_flip()"), ")", collapse = "")))
		} else {
			eval(parse(text = paste("print(qplot(", name.intersect,
				",data=mergedata,geom='histogram')+facet_wrap(~source, scales = 'free_y', ncol = 1)",
				ifelse(is.num, "", "+coord_flip()"), ")", collapse = "")))
		}
		}
		else {
			name.table = matrix(NA, ncol=n, nrow=2)
			for (i in 1:n) {
				name.table[,i] = gt2[[i]][gt4[name.select,1], ]
			}
			name.intersect = as.character(gt4[name.select,2])
			name.class = as.character(gt4[name.select, 3])
			mergedata = data.frame(source = rep(simplifynames(gsub('.csv','',basename(gtfile))),
				rows))
			tmp.num = c()
			for (j in 1:2) {
            for (i in 1:n) {
                if (!is.na(name.table[j,i])) {
                    tmp.num = c(tmp.num, dataset[[i]][, name.table[j,i]])
                } else {
					tmp.num = c(tmp.num, rep(NA, rows[i]))
                }
			}
			mergedata = cbind(mergedata,tmp.num)
			colnames(mergedata)[j+1] = name.intersect[j]
            }
			mergedata = data.frame(mergedata)
			colnames(mergedata)[1]="source"
			mergedata[,2:3] = sapply(mergedata[,2:3],function(avec){as.numeric(as.character(avec))})
			mergedata$source = factor(mergedata$source)
			gbcombo4421 = ggraphics(container = gbcombo442,
				height = 75 * 3 * n,  expand = TRUE)
			if (yscale=="regular y scale") {
				eval(parse(text = paste("print(qplot(", name.intersect[1],",", name.intersect[2],
					",data=mergedata,geom='point',facets=source~.))", sep = "")))
			} else {
				eval(parse(text = paste("print(qplot(", name.intersect[1],",", name.intersect[2],
					",data=mergedata,geom='point')+facet_wrap(~source, scales = 'free', ncol = 1))", sep = "")))
			}
		}
    }

    dict = function(h, ...) {
    #####--------------------------------#####
    ##  dict is the handler of gbcombo432.  ##
	##  (gbutton: Dictionary)               ##
    #####--------------------------------#####
        graphics.off()
			delete(group43, group45)
			group45 <- ggroup(cont=group43,expand = TRUE, use.scrollwindow = TRUE)
			gbcombo44 <- glayout(container = group45,expand = TRUE, use.scrollwindow = TRUE)
        gbcombo44[1, 1, expand = TRUE] = gbcombo443 = gtext(container = gbcombo44, expand = TRUE,
			use.scrollwindow = TRUE)

        name.select = svalue(gt4, index = TRUE)
        if (length(name.select) == 0) {
            gmessage("Please select the variables!")
            gbcombo443[, ] = data.frame(VarName = character(0),
                Level = integer(0), Label = character(0), stringsAsFactors = FALSE)
            return()
        }
        name.table = matrix(nrow = length(name.select), ncol = n)
        for (i in 1:n) {
            name.table[, i] = gt2[[i]][gt4[name.select,1], ]
        }
        name.intersect = as.vector(svalue(gt4))
        name.class = gt4[name.select, 3]

        dictionary = list()
        dictlength = matrix(0, nrow = length(name.intersect),
            ncol = n)
        for (i in 1:length(name.intersect)) {
            dictionary[[i]] = list()
            names(dictionary)[i] = name.intersect[i]
            if (name.class[i] == "factor") {
                for (j in 1:n) {
                  dictionary[[i]][[j]] = levels(factor(dataset[[j]][,
                    name.table[i, j]]))
                  dictlength[i, j] = length(dictionary[[i]][[j]])
                }
            }
        }

        dictlist = list()
        for (i in 1:length(name.intersect)) {
            if (max(dictlength[i, ]) != 0) {
                dictlist[[i]] = matrix(NA, nrow = max(dictlength[i,
                  ]), ncol = n)
                names(dictlist)[i] = name.intersect[i]
                rownames(dictlist[[i]]) = 1:nrow(dictlist[[i]])
				levelintersect = intersect2(dictionary[[i]],dictionary[[i]])
                for (j in 1:n) {
                  dictlist[[i]][1:dictlength[i, j], j] = c(levelintersect$individual[,j],levelintersect$uniq[[j]])
                }
                colnames(dictlist[[i]]) = simplifynames(gsub('.csv','',basename(gtfile)))
            }
            else {
                dictlist[[i]] = "Not a factor"
                names(dictlist)[i] = name.intersect[i]
            }
        }

        if (sum(dictlength) == 0) {
            gmessage("All the variables selected are not factor variables.")
            svalue(gbcombo443) = ""
            return()
        }
        else {
            svalue(gbcombo443) = capture.output(noquote(dictlist))
        }
    }

	change = function(h,...) {
		flagsym = svalue(radio131)

		if (flagsym=="Do not show p-values or flags") {
			newgt4 = gt4[,1:3]
			delete(group42, gt4)
			gt4 <- gtable(newgt4, multiple = T, container = group42, expand = TRUE, chosencol = 2)
			return()
		}

		if (!exists("name_intersection_panel")) {
			mergegui_env$namepanel = nametable
			for (i in 1:n) {
				mergegui_env$namepanel[,i]<-gt2[[i]][,]
			}
			mergegui_env$nameintersection <- gt4[,2]
			mergegui_env$name_intersection_panel <- data.frame(gt4[,],
				scale.rpart(mergegui_env$namepanel,dataset,mergegui_env$nameintersection),
				scale.kstest(mergegui_env$namepanel,dataset,mergegui_env$nameintersection),
				scale.missing(mergegui_env$namepanel,dataset,mergegui_env$nameintersection),
				stringsAsFactors = FALSE)
			colnames(mergegui_env$name_intersection_panel) <- c("Items", "Variables", "Class", "Unit","Dist","Miss")
		} else {
			checknamepanel=c()
			for (i in 1:n) {
				checknamepanel[i]=all(mergegui_env$namepanel[,i]==gt2[[i]][,])
				if (!checknamepanel[i]) mergegui_env$namepanel[,i]<-gt2[[i]][,]
			}
			if (!all(checknamepanel)){
				mergegui_env$nameintersection <- gt4[,2]
				mergegui_env$name_intersection_panel<- data.frame(gt4[,],
					scale.rpart(mergegui_env$namepanel,dataset,mergegui_env$nameintersection),
					scale.kstest(mergegui_env$namepanel,dataset,mergegui_env$nameintersection),
					scale.missing(mergegui_env$namepanel,dataset,mergegui_env$nameintersection),
					stringsAsFactors = FALSE)
			}
		}
		delete(group42, gt4)

		if (flagsym=="Show the flag symbol") {
			flag1 = !is.na(gt4[,4:6])
			flag2 = sapply(gt4[,4:6],function(avec){
				as.numeric(as.character(avec))<0.05
			})
			flag = flag1 & flag2
			newgt4 = mergegui_env$name_intersection_panel
			newgt4[,4]<- as.character(newgt4[,4])
			newgt4[,5]<-as.character(newgt4[,5])
			newgt4[,6]<-as.character(newgt4[,6])
			newgt4[flag[,1],4] <- "X"
			newgt4[flag[,2],5] <- "X"
			newgt4[flag[,3],6] <- "X"
			newgt4[!flag[,1],4] <- ""
			newgt4[!flag[,2],5] <- ""
			newgt4[!flag[,3],6] <- ""
			gt4 <- gtable(newgt4, multiple = T, container = group42,
				expand = TRUE, chosencol = 2)
		} else {
			gt4 <- gtable(mergegui_env$name_intersection_panel, multiple = T,
				container = group42, expand = TRUE, chosencol = 2)
		}

	}

    watchdatafunc = function(h, ...) {
	#####-------------------------------------------------------#####
    ##  watchdatafunc is a function to export the merged dataset.  ##
    ##  For the selected checkboxs, we export the corresponding    ##
	##          variables from all files.                          ##
    ##  The public name for the selected variable is the shortest  ##
	##          name of that variable among different files.       ##
    ##  mergedata is a matrix to save the merged dataset.          ##
    ##  We should write 'xxx.csv'                                  ##
	##          when we export mergedata and save the file.        ##
    #####-------------------------------------------------------#####
		name.select = svalue(gt5, index = TRUE)
        if (length(name.select) == 0) {
            gmessage("Please select the variables!")
            return()
        }
		txtpb = txtProgressBar(min=0,max=1,width = 100,style=3)
        name.table = matrix(nrow = length(name.select), ncol = n)
        for (i in 1:n) {
            name.table[, i] = gt2[[i]][name.select, ]
        }
        name.intersect = as.vector(svalue(gt5))
        name.class = gt5[name.select, 3]
        mergedata = matrix(nrow = sum(rows), ncol = nrow(name.table) +
            1)
        colnames(mergedata) = c("source", name.intersect)
		setTxtProgressBar(txtpb, 0.05)
        for (i in 1:n) {
            tmp = matrix(c(rep(basename(gtfile[i]), rows[i]),
                rep(NA, rows[i] * nrow(name.table))), nrow = rows[i])
            colnames(tmp) = c("source", name.table[, i])
            tmp[, na.omit(name.table[, i])] = as.matrix(dataset[[i]])[,
                na.omit(name.table[, i])]
            mergedata[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i],
                ] = tmp
			setTxtProgressBar(txtpb, (0.05+0.4*i/n))
        }

        mergedatasummary = matrix(c(basename(gtfile), rows),
            nrow = n, dimnames = list(basename(gtfile), c("source",
                "size")))
        for (i in 1:length(name.select)) {
            if (name.class[i] != "NA") {
                if (name.class[i] == "numeric" | name.class[i] ==
                  "integer") {
                  if (name.class[i] == "numeric") {
                    mergedata[, i + 1] = as.numeric(mergedata[,
                      i + 1])
                  }
                  if (name.class[i] == "integer") {
                    mergedata[, i + 1] = as.integer(mergedata[,
                      i + 1])
                  }
                  datasummary = matrix(NA, nrow = n, ncol = 6,
                    dimnames = list(basename(gtfile), paste(name.intersect[i],
                      c("NA#s", "mean", "std", "min", "median",
                        "max"), sep = ".")))
                  for (j in 1:n) {
                    if (!is.na(name.table[i, j])) {
                      tmpdata = dataset[[j]][, name.table[i,
                        j]]
                      datasummary[j, 1] = sum(is.na(tmpdata))
                      datasummary[j, 2] = mean(tmpdata, na.rm = TRUE)
                      datasummary[j, 3] = sd(tmpdata, na.rm = TRUE)
                      datasummary[j, 4] = min(tmpdata, na.rm = TRUE)
                      datasummary[j, 5] = median(tmpdata, na.rm = TRUE)
                      datasummary[j, 6] = max(tmpdata, na.rm = TRUE)
                    }
					#setTxtProgressBar(txtpb, 0.45+(i-1)/n*0.5+j/n*0.5/n)
                  }
                  mergedatasummary = cbind(mergedatasummary,
                    datasummary)
                }
                else {
                  datasummary = matrix(NA, nrow = n, ncol = 9,
                    dimnames = list(basename(gtfile), paste(name.intersect[i],
                      c("NA#s", "levels", "matched_levels", "top1_level",
                        "amount_1", "top2_level", "amount_2",
                        "top3_level", "amount_3"), sep = ".")))
                  matchedlevels = list()
                  for (j in 1:n) {
                    if (!is.na(name.table[i, j])) {
                      if (sum(!is.na(dataset[[j]][, name.table[i,
                        j]])) > 0) {
                        matchedlevels[[j]] = names(table(dataset[[j]][,
                          name.table[i, j]], useNA = "no"))
                      }
                      else {
                        matchedlevels[[j]] = NA
                      }
                    }
                    else {
                      matchedlevels[[j]] = NA
                    }
					#setTxtProgressBar(txtpb, 0.45+(i-1)/n*0.5+j/n*0.25/n)
                  }
                  mtch = intersect2(matchedlevels, matchedlevels)
                  for (j in 1:n) {
                    if (!is.na(name.table[i, j])) {
                      tmpdata = dataset[[j]][, name.table[i,
                        j]]
                      tmptable = sort(table(tmpdata, useNA = "no"),
                        decreasing = TRUE)
                      datasummary[j, 1] = sum(is.na(tmpdata))
                      datasummary[j, 2] = length(tmptable)
                      datasummary[j, 3] = length(mtch$public)
                      datasummary[j, 4] = names(tmptable)[1]
                      datasummary[j, 5] = tmptable[1]
                      if (length(tmptable) > 1) {
                        datasummary[j, 6] = names(tmptable)[2]
                        datasummary[j, 7] = tmptable[2]
                      }
                      if (length(tmptable) > 2) {
                        datasummary[j, 8] = names(tmptable)[3]
                        datasummary[j, 9] = tmptable[3]
                      }
                    }
					#setTxtProgressBar(txtpb, 0.45+(i-0.5)/n*0.5+j/n*0.25/n)
                  }
                  mergedatasummary = cbind(mergedatasummary,
                    datasummary)
                }
            }
            else {
                mergedatasummary = cbind(mergedatasummary, matrix(NA,
                  ncol = 1, nrow = n, dimnames = list(NULL, name.intersect[i])))
            }
			setTxtProgressBar(txtpb, 0.45+i/n*0.5)
        }

        setTxtProgressBar(txtpb, 1)
        if (!is.na(gf <- gfile(type = "save"))) {
			if (regexpr("\\.csv$",gf) %in% c(-1,1)) {
				gf = paste(gf,".csv",sep="")
			}
            write.csv(mergedata, file = gf, row.names = F)
            summarylocation = sub("\\.csv$", "_summary.csv", gf)
            write.csv(mergedatasummary, file = summarylocation,
                row.names = F)
            gmessage("The files are merged!")
        }
    }


    #####-------------------------------------------------------------------#####
    ##  Import the selected files.                                             ##
    ##  'dataset' is a list to save all the data from different files.         ##
    ##  'rows' is  a vector to save the number of observations for each file.  ##
    ##  'vname' is  a list to save the original colnames of the dataset.       ##
    ##  'simplifiedname' is  a list to save the simplified name.               ##
	##          (delete the filenames in the colnames if they have)            ##
    ##  'vname' & 'simplifiedname' are 1-1 projections,                        ##
	##          although 'simplifiedname' haves repeated names.                ##
    #####-------------------------------------------------------------------#####
    dataset <- list()
    simplifiedname <- list()
    vname <- list()
	# gt = c("E:\\study\\2010 Fall\\research\\Novartis\\health\\brfss09metric.csv",
	#	"E:\\study\\2010 Fall\\research\\Novartis\\health\\brfss09std.csv")
	# gt = c("E:\\study\\2010 Fall\\research\\Novartis\\oil spill\\Brooks-McCall Cruise01 05-08-2010_QC 2010132.csv",
	#   "E:\\study\\2010 Fall\\research\\Novartis\\oil spill\\Brooks-McCall Cruise02 05-15-2010_QC 2010138.csv")

    if (length(svalue(gt)) == 0) {
        n <- length(gt[])
        gtfile <- gt[]
    } else {
        n <- length(svalue(gt))
        gtfile <- svalue(gt)
    }

	rows <- rep(0, n)
	if (n<2) {
	    warning('The input data set is not enough. More files are needed.')
		return()
	}

    for (i in 1:n) {
        dataset[[i]] <- read.csv(file = gtfile[i], head = TRUE)
        rows[i] <- nrow(dataset[[i]])
        vname[[i]] <- colnames(dataset[[i]])
        J = gregexpr(".csv.", vname[[i]])
        loc = c()
        for (j in 1:length(vname[[i]])) {
            loc[j] = max(J[[j]]) + 5
        }
        loc[which(loc == 4)] = 1
        simplifiedname[[i]] = substring(vname[[i]], loc)
    }

    #####----------------------------------------------------#####
    ##  Now we are going to generate two stuffs:                ##
	##          'nameintersect' and 'nametable'.                ##
    ##  'nametable' is a matrix with all matched and            ##
	##          unmatched variable names ('vname').             ##
    ##  'nameintersect' is a vector of mutually exclusive       ##
	##          and collectively exhaustive names,              ##
	##          with 'simplifiedname' at the beginning,         ##
	##          and the special 'vname' at the end.             ##
	##  length(nameintersect) == nrow(nametable)                ##
    ##  'Part1-1-' is the intersection for all the n files.     ##
    ##  'Part(n-i+1)' is the intersection for all combination   ##
    ##	        of the n-i+1 files.                             ##
    ##  'Partn-i-' is the left part of the i-th files,          ##
	##          it cannot be intersected with any other files.  ##
    #####----------------------------------------------------#####
    a = intersect2(vname, simplifiedname)
    nameintersect <- a$public
    tmpuniq = a$uniq
    tmpsimpleuniq = a$simpleuniq
    nametable <- a$individual
    if (nrow(nametable) != 0) {
        rownames(nametable) <- paste("Part1-1-", 1:nrow(nametable),
            sep = "")
    }

    for (i in max((n - 1), 2):2) {
        combnmatrix = combn(1:n, i)
        for (j in 1:ncol(combnmatrix)) {
            tmpintersect = intersect2(tmpuniq[combnmatrix[, j]],
                tmpsimpleuniq[combnmatrix[, j]])
            tmptable = matrix(NA, ncol = n, nrow = length(tmpintersect$public))
            tmptable[, combnmatrix[, j]] = tmpintersect$individual
            if (nrow(tmptable) != 0) {
                rownames(tmptable) = paste("Part", n - i + 1,
                  "-", j, "-", 1:length(tmpintersect$public),
                  sep = "")
            }
            nametable <- rbind(nametable, tmptable)
            nameintersect <- c(nameintersect, tmpintersect$public)
            tmpuniq[combnmatrix[, j]] = tmpintersect$uniq
            tmpsimpleuniq[combnmatrix[, j]] = tmpintersect$simpleuniq
        }
    }
    nameintersect <- c(nameintersect, unlist(tmpuniq))

    for (i in 1:n) {
        tmptable = matrix(NA, ncol = n, nrow = length(tmpuniq[[i]]))
        tmptable[, i] = tmpuniq[[i]]
        if (nrow(tmptable) != 0) {
            rownames(tmptable) = paste("Part", n, "-", i, "-",
                1:nrow(tmptable), sep = "")
        }
        nametable <- rbind(nametable, tmptable)
    }
    colnames(nametable) <- simplifynames(gsub('.csv','',basename(gtfile)))

	#####-------------------------------#####
    ##  New window for matching variables  ##
    #####-------------------------------#####
 	combo2 = gwindow("Matched Variables", visible = T, width = 900,
        height = 600)
    tab = gnotebook(container = combo2)

	#####----------------------------------------------#####
    ##  In the first tab we can:                          ##
    ##  (1) . ##
    ##  (2) .    ##
    #####----------------------------------------------#####
	group11 = ggroup(horizontal = FALSE, cont = tab, label = "Preferences",
        expand = T)
    frame12 = gframe("Scaling of histograms",container = group11, horizontal = FALSE)
	radio121 <- gradio(c("regular y scale","relative y scale"),container = frame12)
	frame13 = gframe("Flag for variables",container = group11, horizontal = FALSE)
	radio131 <- gradio(c("Do not show p-values or flags","Show p-values","Show the flag symbol"), container = frame13, handler=change)

    #####----------------------------------------------#####
    ##  In the first tab we can:                          ##
    ##  (1) Switch the variable names in the same gtable. ##
    ##  (2) Go back or go forth or reset the matching.    ##
    #####----------------------------------------------#####
    group21 = ggroup(horizontal = FALSE, cont = tab, label = "Matching",
        expand = T)
    group22 = ggroup(container = group21, use.scrollwindow = TRUE,
        expand = T)
    group2 = list()
    gt2 <- list()
    mergegui_env$hstry1 <- list()
    mergegui_env$hstry1[[1]] <- nametable
    mergegui_env$hstry2 <- list()
    mergegui_env$hstry2[[1]] <- nameintersect
    mergegui_env$hstry3 <- list()
    mergegui_env$hstry3[[1]] <- var.class(nametable,dataset)
    mergegui_env$idx <- 1
    mergegui_env$redo.indicate <- 0
    for (i in 1:n) {
        group2[[i]] = ggroup(horizontal = FALSE, container = group22,
            expand = T)
        gt2[[i]] <- gtable(nametable[, i, drop = F], container = group2[[i]],
            expand = TRUE)
        tag(gt2[[i]], "prev.idx") <- svalue(gt2[[i]], index = TRUE)
        tag(gt2[[i]], "toggle") <- FALSE
        tag(gt2[[i]], "idx") <- i
        addhandlerclicked(gt2[[i]], handler = function(h, ...) {
            gt.tmp = h$obj
            prev.idx = tag(gt.tmp, "prev.idx")
            if (length(prev.idx) == 1 && tag(gt.tmp, "toggle")) {
                tmp = gt.tmp[prev.idx, ]
                gt.tmp[prev.idx, ] = svalue(gt.tmp)
                gt.tmp[svalue(gt.tmp, index = TRUE), ] = tmp
                mergegui_env$idx <- mergegui_env$idx + 1
                mergegui_env$hstry1[[mergegui_env$idx]] <- mergegui_env$hstry1[[mergegui_env$idx - 1]]
                mergegui_env$hstry1[[mergegui_env$idx]][, tag(gt.tmp, "idx")] <- gt.tmp[]
                if (tag(gt.tmp, "idx") == 1) {
                  tmpgt4 = gt4[prev.idx, 2:3]
                  gt4[prev.idx, 2:3] = gt4[svalue(gt.tmp, index = TRUE), 2:3]
                  gt4[svalue(gt.tmp, index = TRUE), 2:3] = tmpgt4
                }
                mergegui_env$hstry2[[mergegui_env$idx]] <- gt4[, 2]
                mergegui_env$hstry3[[mergegui_env$idx]] <- gt4[, 3]
                gt5[, 2] <- gt4[, 2]
                gt5[, 3] <- gt4[, 3]
                mergegui_env$redo.indicate <- 0
            }
            tag(gt.tmp, "toggle") = !tag(gt.tmp, "toggle")
            tag(gt.tmp, "prev.idx") = svalue(gt.tmp, index = TRUE)
        })
    }
    group23 <- ggroup(container = group21)
    gbcombo21 <- gbutton("Undo", container = group23, handler = undo,
        expand = TRUE)
    gbcombo22 <- gbutton("Redo", container = group23, handler = redo,
        expand = TRUE)
    gbcombo23 <- gbutton("Reset", container = group23, handler = reset,
        expand = TRUE)

    #####------------------------------------------------#####
    ##  In the second tab we can:                            ##
    ##  (1) Watch and change the name or type of variables. ##
    ##  (2) Numeric or graphic summary.                     ##
    ##  (3) Dictionary for factor variables.                ##
    #####------------------------------------------------#####
    group41 = ggroup(cont = tab, label = "Summary", expand = T)
    group42 <- ggroup(container = group41, use.scrollwindow = TRUE,
        expand = T)
    name.inter <- data.frame(1:length(nameintersect), nameintersect,
        var.class(nametable,dataset), stringsAsFactors = FALSE)
    colnames(name.inter) = c("Items", "Variables", "Class")
    gt4 <- gtable(name.inter, multiple = T, container = group42,
        expand = TRUE, chosencol = 2)
    addhandlerdoubleclick(gt4, handler = VariableOptions)
    group43 <- ggroup(horizontal = FALSE, container = group41,
        expand = TRUE)
    group44 = ggroup(horizontal = TRUE, container = group43)
    gbcombo431 <- gbutton("Numeric summary", container = group44,
        handler = smmry, expand = TRUE)
    gbcombo432 <- gbutton("Graphical summary", container = group44,
        handler = graph, expand = TRUE)
    gbcombo433 <- gbutton("Dictionary", container = group44, handler = dict,
        expand = TRUE)
	group45 <- ggroup(container = group43, expand = TRUE, use.scrollwindow = TRUE)

    #####------------------------------------------------#####
    ##  In the third tab we can:                           ##
    ##  (1) Select all or none variables.                   ##
    ##  (2) Export the data.                                ##
    #####------------------------------------------------#####
    group51 = ggroup(container = tab, label = "Export", expand = T)
    group52 = ggroup(container = group51, use.scrollwindow = TRUE,
        expand = T)
	Matched = substr(rownames(nametable),5,regexpr('-',rownames(nametable))-1)
	FileMatched = as.character((n+1)-as.integer(Matched))
    gt5 <- gtable(data.frame(name.inter[,1:3],FileMatched), multiple = T, container = group52,
        expand = TRUE, chosencol = 2)
	addhandlerclicked(gt5,handler=function(h,...){
		svalue(gbcombo57) = paste("Currently you select",length(svalue(gt5)),"variables.",sep=" ")
	})
    group53 = ggroup(horizontal = FALSE, container = group51,
        expand = TRUE)
    gbcombo51 <- gbutton("Select All", container = group53, handler = function(h,
        ...) {
        svalue(gt5, index = TRUE) = 1:length(nameintersect)
		svalue(gbcombo57) = paste("Currently you select all",length(nameintersect),"variables.",sep=" ")
        focus(gt5)
    })
    gbcombo52 <- gbutton("Clear All", container = group53, handler = function(h,
        ...) {
        svalue(gt5) = NULL
		svalue(gbcombo57) = "Currently you select 0 variable."
    })
    gbcombo55 <- gbutton("Export the matched data", container = group53,
        handler = watchdatafunc)
	gbcombo56 <- glabel(paste("The complete merged data have ",sum(rows)," rows and ",
		length(nameintersect)," columns."),cont=group53)
	gbcombo57 <- glabel(paste("Currently you select 0 variable."),cont=group53)

	svalue(tab)=1
}

mergeID = function(h, ...) {
	####################################################
	# mergeID is a function to merge the observations. #
	####################################################

	watchIDfunc = function(h, ...) {
    #####------------------------------------------------------------------------------------#####
    ##  watchIDfunc is a function to export the merged dataset.                                 ##
    ##  key is a vector of the selected primary keys. We checked the validity of the key first. ##
    ##  keyID is the merged ID for the observations.                                            ##
    ##  mergeIDdata is the matrix that merged all the files by the keyID.                       ##
    ##  We should write 'xxx.csv' when we export mergeIDdata and save the file.                 ##
    #####------------------------------------------------------------------------------------#####
        keyID = c()
        vcolumn = rep(0, n)
        key = c()
        for (i in 1:n) {
            key[i] = svalue(gt3[[i]])
            if (sum(duplicated(as.character(dataset[[i]][, key[i]]))) >
                0) {
                gmessage(paste(key[i], "could not be the primary key for",
                  basename(gtfile[i]), "because it has repeated items. Please choose another key."))
                return()
            }
            keyID = union(keyID, dataset[[i]][, key[i]])
            vcolumn[i] = length(vname[[i]]) - 1
        }
        mergeIDdata = matrix(NA, nrow = length(keyID), ncol = sum(vcolumn) +
            1, dimnames = list(keyID))
        mergeIDdata[, 1] = keyID
        mergeIDcolnames = c(key[1], 1:sum(vcolumn) + 1)
        for (i in 1:n) {
            mergeIDdata[as.character(dataset[[i]][, key[i]]),
                cumsum(c(2, vcolumn))[i]:cumsum(c(1, vcolumn))[i +
                  1]] = as.matrix(dataset[[i]][, setdiff(vname[[i]],
                key[i])])
            mergeIDcolnames[cumsum(c(2, vcolumn))[i]:cumsum(c(1,
                vcolumn))[i + 1]] = paste(basename(gtfile[i]),
                ".", colnames(dataset[[i]][, setdiff(vname[[i]],
                  key[i])]), sep = "")
        }
        colnames(mergeIDdata) = mergeIDcolnames

        if (!is.na(gf <- gfile(type = "save"))) {
			if (regexpr("\\.csv$",gf) %in% c(-1,1)) {
				write.csv(mergeIDdata, file = paste(gf,".csv",sep=""), row.names = FALSE)
			} else {
				write.csv(mergeIDdata, file = gf, row.names = FALSE)
			}
			gmessage("The files are merged!")
		}
    }

    #####-----------------------------------------------------------#####
    ##  Input the selected files.                                      ##
    ##  dataset is a list to save all the data from different files.   ##
    ##  vname is  a list to save the original colnames of the dataset. ##
    #####-----------------------------------------------------------#####
    dataset <- list()
    vname <- list()
    simplifiedname <- list()
    if (length(svalue(gt)) == 0) {
        n <- length(gt[])
        gtfile <- gt[]
    }
    else {
        n <- length(svalue(gt))
        gtfile <- svalue(gt)
    }
    rows <- rep(0, n)
    for (i in 1:n) {
        dataset[[i]] <- read.csv(file = gtfile[i], head = T)
        rows[i] <- nrow(dataset[[i]])
        vname[[i]] <- colnames(dataset[[i]])
        J = gregexpr(".csv.", vname[[i]])
        loc = c()
        for (j in 1:length(vname[[i]])) {
            loc[j] = max(J[[j]]) + 5
        }
        loc[which(loc == 4)] = 1
        simplifiedname[[i]] = substring(vname[[i]], loc)
    }

    a = intersect2(vname, simplifiedname)
    tmpuniq = a$uniq
    tmpnametable = a$individual
    tmpvname = list()
    for (i in 1:n) {
        tmpvname[[i]] = c(tmpnametable[, i], tmpuniq[[i]])
    }

    #####-----------------------------------------#####
    ##  In this GUI we can:                          ##
    ##  Select the primary keys for different files. ##
    #####-----------------------------------------#####
    combo3 <- gwindow("Matched Primary Key", visible = TRUE)
    group3 <- ggroup(horizontal = FALSE, container = combo3)
    gt3 <- list()

    for (i in 1:n) {
        gl3 = glabel(paste("Select the Primary Key from", basename(gtfile[i])),
            container = group3)
        gt3[[i]] <- gcombobox(tmpvname[[i]], container = group3,
            expand = T)
    }
    gbcombo3 = gbutton("Match by the Key", container = group3,
        handler = watchIDfunc)

}

	#####---------------------#####
	##  First GUI:  Open files.  ##
	#####---------------------#####

	combo <- gwindow("Combination", visible = TRUE)
	group <- ggroup(horizontal = FALSE, container = combo)
	f.list <- matrix(nrow = 0, ncol = 1, dimnames = list(NULL,
		"File"))

	gt <- gtable(f.list, multiple = T, container = group,
		expand = T)
	gb1 <- gbutton("Open", container = group, handler = function(h,
		...) gt[, ] = na.omit(rbind(gt[, , drop = FALSE], matrix(if (.Platform$OS.type != 'windows')
                                                   file.choose() else choose.files(), dimnames = list(NULL,
		"File")))))
	gb2 <- gbutton("Match the Variables", container = group,
		handler = mergefunc)
	gb3 <- gbutton("Match by the Key", container = group,
		handler = mergeID)
}
