\name{MergeGUI}
\alias{MergeGUI}
\title{The Merging GUI.
This function will start with an starting interface, allowing 1)
selecting several data files; 2) doing the next command with more
than one files. There are two commands which could be selected:
match the variables, match the cases by the key variable. In the
matching-variable interface the user can 1) check the matching of
the variables among files and switch the variable names if they
are wrongly matched; 2) look at the numerical and graphical
summaries for the selected variables, or the dictionary for
selected factor varibles; 3) observe the misclassification rate,
KS-test p-values and Chi-square test p-values for each variable,
which helps to determine whether any transformation is needed for
the variable; (For each variable, the user may want to know
whether it could distinguish the sources correctly. So the
misclassification rate is calculated through the tree model.
KS-test is used to check whether any variable has different
distributions for different sources. And the Chi-square test is
useful when the user is interested in the pattern of missing
values among the sources.) 4) change the name or class for any
variable; 5) export the merged dataset and the summary for it. In
the matching-case interface the user can determine a primary key
for each data file and then merge the cases by the key.}
\usage{
MergeGUI(..., filenames = NULL, unit = TRUE, distn = TRUE, miss = TRUE)
}
\arguments{
  \item{...}{names of the data frames to read}

  \item{filenames}{A vector of csv file names with their
  full paths.}

  \item{unit}{whether the test of the difference among the
  group centers is on or off}

  \item{distn}{whether the test of the difference among the
  group distributions is on or off}

  \item{miss}{whether the test of the difference among the
  group missing patterns is on or off}
}
\value{
  NULL
}
\description{
  The merging GUI consists of four tabs. In the preferences
  tab, user can choose whether the numerical p-values or
  the flag symbols are displayed in the summary tab;
  whether the y-scales are free for different data files
  when drawing the plots faceted by the sources. In the
  checking tab, each data file has a list of variable
  names, and the GUI will automatically arrange the order
  of variable names to align the same names in one row. The
  user can switch the order of the variables in one file's
  list. It is possible to undo, redo, or reset the
  matching. In the summary tab, there is a list of variable
  names on the left which corresponds to the checking tab.
  The misclassification rate, KS-test p-values and
  Chi-square test p-values for each variable may also be
  presented with the variable names. On the top right there
  are three buttons: Numeric summary, Graphical summary,
  and Dictionary. And the results could be shown below the
  buttons. For the graphical summary, histogram or barchart
  will be shown if a single variable is selected. A
  scatterplot will be drawn if two numeric or two factor
  varaibles are chosen. Side-by-side boxplots will be
  presented when one numeric and one factor varaibles are
  selected. A parallel coordinate plot is shown when all
  the variables selected are numeric and there are more
  than two variables. If more than two variables are chosen
  but the classes of the variables are mixed, i.e. some are
  numeric, some are factor or character, then histograms
  and barcharts will be drawn individually. All the plots
  are facetted by the source. In the export tab the user
  could select all or none variables by click the buttons
  or choose several varaibles by Ctrl+Click. Then the
  export button will export the merged data and the numeric
  summaries of the selected variables into two csv files.
}
\examples{
if (interactive()) {
    MergeGUI()
    
    csvnames = list.files(system.file("doc", package = "MergeGUI"), 
        pattern = "\\\\.csv$")
    files = system.file("doc", csvnames, package = "MergeGUI")
    MergeGUI(filenames = files)
    
    data(iris)
    setosa = iris[iris$Species == "setosa", 1:4]
    versicolor = iris[iris$Species == "versicolor", 1:4]
    virginica = iris[iris$Species == "virginica", 1:4]
    MergeGUI(setosa, versicolor, virginica)
}
}
\author{
  Xiaoyue Cheng <\email{xycheng@iastate.edu}>
}
