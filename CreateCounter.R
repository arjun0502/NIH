#Arjun Jain 
#Helper function to create a counter and increment
#Used in CalculateStainingIndex, FindCounters, and AdvancedCalculateStainingIndex 

CreateCounter <- function(curr.count) {
  list(
    increment = function(amount) {
      curr.count <<- curr.count + amount
    },
    value = function() {
      return(curr.count)
    }
  )
}