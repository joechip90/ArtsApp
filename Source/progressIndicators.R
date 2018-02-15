## 1. ------ DEFINE A FUNCTION TO CREATE A PROGRESS BAR UPDATER FUNCTION ------
#' @title Create a Function to Update a Progress Bar
#' 
#' @description 
#' This function creates another function that can be used to update a progress
#' bar.  Useful to pass to functions that may take a long time to run and needs
#' to give user feedback on the progress.
#' 
#' @param progressOb An object representing a progress indicators.  Either this
#' can be the output from the \code{\link[utils::txtProgressBar]{txtProgressBar}}
#' function or from the \code{\link[shiny::Progress]{Progress}} function.
#' 
#' @return A function that can be called to update the progress object given in
#' \code{progressOb}.  The function can take the following arguments:
#' \describe{
#'     \item{\code{message}}{A message to display to the user about the currently
#'     executing process}
#'     \item{\code{details}}{A message to display to the user about the details
#'     of the currently executing process}
#'     \item{\code{value}}{A value between the minimum and maximum values set in
#'     the progress object}
#' }
#' 
#' @examples
#' # Initialise a progress bar to use
#' examProgress <- txtProgressBar()
#' # Create a function to update the progress bar
#' examProgressUpdate <- createProgressUpdater(examProgress)
#' # Use the function to update the progress bar
#' examProgressUpdate(
#'     message = "This is a test progress update",
#'     details = "Lets say it is 50% done",
#'     value = 0.5)
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[utils::txtProgressBar]{txtProgressBar}}
#' \code{\link[shiny::Progress]{Progress}}
#' @export
#' 
createProgressUpdater <- function(progressOb) {
	### 1.1 ==== Sanity test the inputs ====
	# Determine the type of progress bar initialised
	if(class(progressOb) == "txtProgressBar") {
		inProgressType <- "utils"
	} else if(all(class(progressOb) == c("Progress", "R6"))) {
		inProgressType <- "shiny"
	}
	outFunc <- NULL
	### 1.2 ==== Create an updater function for shiny progress bars ====
	if(as.character(inProgressType) == "shiny") {
		outFunc <- function(message, details, value, valProgressOb = progressOb) {
			## 1.2.1 Sanity test the inputs ----
			# Retrieve the updater message
			inMessage <- message
			if(!is.null(inMessage)) {
				# Convert the message to a character vector
				inMessage <- tryCatch(as.character(inMessage), error = function(err) {
					c("")
				})
				# Swap out any NAs to ""
				inMessage[is.na(inMessage)] <- ""
				# Collapse the message into one element
				inMessage <- paste(inMessage, collapse = " ")
			}
			# Retrieve the updater details
			inDetails <- details
			if(!is.null(inDetails)) {
				# Convert the details to a character vector
				inDetails <- tryCatch(as.character(inDetails), error = function(err) {
					c("")
				})
				# Swap out any NAs to ""
				inDetails[is.na(inDetails)] <- ""
				# Collapse the message into one element
				inDetails <- paste(inDetails, collapse = " ")
			}
			# Retrieve the updater value
			inValue <- value
			if(!is.null(inValue)) {
				# Convert the value to a numeric vector
				inValue <- tryCatch(as.numeric(inValue), error = function(err) {
					c()
				})
				# Set the input value to NULL if the input is NA, NaN, non-finite or length 0
				if(length(inValue) <= 0 || is.na(inValue[1]) || !is.finite(inValue[1])) {
					inValue <- NULL
				} else {
					inValue <- inValue[1]
				}
			}
			## 1.2.2 Update the progress object ----
			# Update the message if requested
			tryCatch(valProgressOb$set(message = inMessage, detail = inDetails, value = inValue), error = function(err) {
				stop(paste("unable to update progress indicator:", err, sep = " "))	
			})
		}
	### 1.3 ==== Create an updater function for utils progress bars ====	
	} else if(as.character(inProgressType) == "utils") {
		outFunc <- function(message, details, value, valProgressOb = progressOb) {
			## 1.3.1 Sanity test the inputs ----
			# Retrieve the updater message
			inMessage <- message
			if(!is.null(inMessage)) {
				# Convert the message to a character vector
				inMessage <- tryCatch(as.character(inMessage), error = function(err) {
					c("")
				})
				# Swap out any NAs to ""
				inMessage[is.na(inMessage)] <- ""
				# Collapse the message into one element
				inMessage <- paste(inMessage, collapse = " ")
			}
			# Retrieve the updater details
			inDetails <- details
			if(!is.null(inDetails)) {
				# Convert the details to a character vector
				inDetails <- tryCatch(as.character(inDetails), error = function(err) {
					c("")
				})
				# Swap out any NAs to ""
				inDetails[is.na(inDetails)] <- ""
				# Collapse the message into one element
				inDetails <- paste(inDetails, collapse = " ")
			}
			# Retrieve the updater value
			inValue <- value
			if(!is.null(inValue)) {
				# Convert the value to a numeric vector
				inValue <- tryCatch(as.numeric(inValue), error = function(err) {
					c()
				})
				# Set the input value to NULL if the input is NA, NaN, non-finite or length 0
				if(length(inValue) <= 0 || is.na(inValue[1]) || !is.finite(inValue[1])) {
					inValue <- NULL
				} else {
					inValue <- inValue[1]
				}
			}
			## 1.3.2 Update the progress object ----
			if(is.null(inValue)) {
				inValue <- NA
			}
			tryCatch(setTxtProgressBar(valProgressOb, value = inValue, title = inMessage, label = inDetails), error = function(err) {
				stop(paste("unable to update progress indicator:", err, sep = " "))
			})
		}
	} else {
		stop("unknown progress indicator object type")
	}
	outFunc
}

## 2. ------ DEFINE A FUNCTION TO UPDATE THE PROGRESS BAR ------
#' @title Update Progress Bar
#' 
#' @description 
#' This function calls a progress update function and passes it the inputs
#' provided as parameters.  This functions is largely only used as a helper
#' functions inside functions that take the output from
#' \code{\link[createProgressUpdater]{createProgressUpdater}} as an argument-
#' 
#' @param updateFunc A function created using
#' \code{\link[createProgressUpdater]{createProgressUpdater}}.
#' @param value Numeric scalar.  A value between the minimum and maximum
#' values set in the progress object.
#' @param details Character scalar. A message to display to the user about
#' the details of the currently executing process.
#' @param message Character scalar. A message to display to the user about
#' the currently executing process.
#' 
#' @return For progress bars created by \code{\link[utils::txtProgressBar]{txtProgressBar}},
#' the function returns the output from \code{\link[utils::txtProgressBar]{setProgressBar}}.
#' For progress bars created by \code{\link[shiny::Progress]{Progress}}, the
#' function returns the output from the \code{set} method of the respective
#' \code{\link[shiny::Progress]{Progress}} object.
#' 
#' @examples
#' # Initialise a progress bar to use
#' examProgress <- txtProgressBar()
#' # Create a function to update the progress bar
#' examProgressUpdate <- createProgressUpdater(examProgress)
#' # Use the function to update the progress bar
#' updateProgressBar(
#'     updateFunc = examProgressUpdate,
#'     message = "This is a test progress update",
#'     details = "Lets say it is 50% done",
#'     value = 0.5)
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[utils::txtProgressBar]{txtProgressBar}}
#' \code{\link[shiny::Progress]{Progress}}
#' @export
#' 
updateProgressBar <- function(updateFunc, value, details = NULL, message = NULL) {
	progressOut <- NULL
	if(!is.null(updateFunc)) {
		# Coerce the updateFunc input to a function
		inUpdateFunc <- tryCatch(as.function(updateFunc), error = function(err) {
			stop(paste("invalid progress bar update function:", err, sep = " "))
		})
		# Use the updateFunc to update the progress bar
		progressOut <- tryCatch(inUpdateFunc(value = value, details = details, message = message), error = function(err) {
			stop(paste("error encountered during application of progress bar update function:", err, sep = " "))
		})
	}
	progressOut
}