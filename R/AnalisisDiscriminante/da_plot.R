da_ggplot <- function(model, data, class, 
                                predict_type = "class",
                                resolution = 1000) {

  data <- data[,1:2]
  cn <- colnames(data)
  
  k <- length(unique(class))
  data$Col <-class
  
  gg <- ggplot(data=data, aes(data[,1], data[,2],color=Col))+ 
    geom_point(show.legend = F)+
    theme_minimal()+
    labs(x='',y='')
  
  # make grid
  r <- sapply(data[, 1:2], range, na.rm = TRUE)
  xs <- seq(r[1, 1], r[2, 1], length.out = resolution)
  ys <- seq(r[1, 2], r[2, 2], length.out = resolution)
  
  g <- cbind(rep(xs, each = resolution), rep(ys, time = resolution))
  
  colnames(g) <- colnames(r)
  
  g <- as.data.frame(g)
  
  # Predict 
  p <- predict(model, g, type = predict_type)
  
  if(is.list(p)){
    p <- p$class
  } 
  g$col<- as.integer(as.factor(p))
  
  gg + geom_contour(aes(x = g[,1], y = g[,2], z = g[,3]), data = g, 
                    col="black",linewidth=.25)
}
