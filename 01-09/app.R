
library (shiny)
library(shinythemes)
library(ggplot2)
library(SciViews)
library(Ryacas)
Sys.setlocale(category = "LC_ALL", locale = "Spanish")

##FUNCIONES
feval <- function(f, ...) do.call(f, list(...))

#Integración
Rtrapecio <- function(F, a,b) {
  h<-(b-a)/1
  
  result <- (h/2)*(feval(F, a)+feval(F, b))
    
  return(round(result,3))
}

Rsimpson <- function(F, a, b) {
  h <-(b-a)/2
  
  result <- (h/3)*(feval(F, a)+4*feval(F, a+h)+feval(F, b))
  
  return(round(result,3))
}

Rsimpson3oct <- function(F, a, b) {
  h <-(b-a)/3
  
  result <- ((3*h)/8)*(feval(F, a)+3*feval(F, a+h)+3*feval(F, a+2*h)+feval(F, b))
  
  return(round(result,3) )
}



RPuntoMedio <-function(F, a, b) {
  h <-(b-a)/2
  x0 <-a+h
  result <- 2*h*feval(F,x0)
    
  return (round(result,5))
}

fabierta2P <-function(F, a, b) {
  h <-(b-a)/3
  x0 <-a+h
  x1 <-b-h
  result <-((3*h)/2)*(feval(F,x0)+feval(F,x1))
  
  return (round(result,5))
}

fabierta3P <-function(F, a, b) {
  h <-(b-a)/4
  x0 <-a+h
  x1 <-x0+h
  x2 <-b-h
  result <-((4*h)/3)*(2*feval(F,x0)-feval(F,x1)+2*feval(F,x2))
  
  return (round(result,5))
}


RtrapecioC <-function(F, a, b,N) {
  h <-(b-a)/N
  
  result <-feval(F, a)
  sum<-0
  for(j in 2: N-1){
    sum<- sum+ feval(F, a+j*h)
  }
  sum<-2*sum
  result<- (h/2)*(result+sum+ feval(F, b))
  
  return (result)
}


RsimpsonC <-function(F, a, b,N) {
  h <-(b-a)/N
  if(a>=0){
    result <-feval(F, a)
    sum<-0
      for(j in 2: (N/2)-1){
        sum<- sum+ feval(F, a+2*j*h)
      }
    
    sum<-2*sum
    sum2<-0
    
      for(j in 2: (N/2)){
        sum2<- sum2+ feval(F, a+((2*j)-1)*h)
      }
    
    sum2<-4*sum2
    result<- (h/3)*(result+sum+sum2+ feval(F, b))
  }
  else{
    result <-feval(F, a)
  
    sum<-0
    
    for(j in 2: (N/2)-1){
    
      sum<- sum+ feval(F, a+2*j*h)
      
    }
    
    sum<-2*sum
    sum2<-0
    for(j in 1: (N/2)){
      
      sum2<- sum2+ feval(F, a+((2*j)-1)*h)
     
    }
    
    sum2<-4*sum2
    result<- (h/3)*(result+sum+sum2+ feval(F, b))
  }
  return (result)
}


RPuntoMedioC <-function(F, a, b,N) {
  h <-(b-a)/(N+2)

  sum<-0
  if(a<0){
    for(j in  (N/2):0){

      sum<- sum+ feval(F, 2*j*(a/N))

    }
  }else{
  for(j in 0: (N/2)){
  
    sum<- sum+ feval(F, 2*j*(b/N))

  }
  }
 
  result<- 2*h*(sum)

  return (result)
}



#Derivación
derivC <- function(F, x0, h) {
  
  d <- (feval(F, x0 + h) - feval(F, x0 - h))/(2 * h)
  
  return(round(d,3))
  
}
derivI <- function(F, x0, h) {
  
  d <- (feval(F, x0) - feval(F, x0 - h))/( h)
  
  return(round(d,3))
  
}

derivD <- function(F, x0, h) {
  
  d <- (feval(F, x0 + h) - feval(F, x0))/( h)
  
  return(round(d,3))
  
}

#Interpolación

## hermite
# p_deriv <- polyder(p_real)
# dy <- polyval(p_deriv, x)
# p_hermite <- hermite(x, y, dy)
# sum(abs(polyval(p_hermite, x) - y))
# sum(abs(polyval(p_deriv, x) - dy))

#' @export
#' @importFrom pracma zeros
hermite_table <- function(x, y, dy, decimals = 0) {
  
  n <- length(x)
  
  xx <- c(x, x)
  
  if (decimals > 0)
    xx <- xx %>% round(digits = decimals)
  
  yy <- c(y, y)
  if (decimals > 0)
    yy <- yy %>% round(digits = decimals)
  
  if (decimals > 0)
    dy <- dy %>% round(digits = decimals)
  
  
  ii <- order(xx)
  xx <- xx[ii]
  yy <- yy[ii]
  
  z <- pracma::zeros(2 * n, 2 * n)
  
  z[, 1] <- t(yy)
  
  
  for (col in 2:(2 * n)) {
    
    for (fil in 1:(2 * n - col + 1)) {
      
      if (xx[fil + col - 1] != xx[fil]) {
        
        z[fil, col] <- (z[fil + 1, col - 1] - z[fil,
                                                col - 1])/(xx[fil + col - 1] - xx[fil])
        
        if (decimals > 0)
          z[fil, col] <- z[fil, col] %>% round(digits = decimals)
        
        
      }
      
    }
    
    if (col == 2) {
      
      z[seq(1, 2 * n, by = 2), 2] <- t(dy)
      
      
    }
    
  }
  
  return(z)
  
}

polinomioInt<- function(fila1,columna1){
  result<-""
 
  for(fil in 2:dim(fila1)[2]){
  
    if(fil==2){
      result<-paste(result,round(fila1[fil],3), sep = " ")
    }else if(fil!=1) {
      result<-paste(result,"+(",round(fila1[fil],3),sep = " ")
     for(numF in 3:fil){
       if(columna1[numF-2]<0){
         result<-paste(result,"(x +",round(abs(columna1[numF-2]),3),")", sep = " ")
       }else{
        result<-paste(result,"(x -",round(columna1[numF-2],3),")", sep = " ")
       }
     }
     result<-paste(result,")",sep = " ")
    }
    
   
  }
  return (result)
}
polinomioInt2<- function(fila1,columna1,x){
  result<-""
  
  for(fil in 2:dim(fila1)[2]){
    
    if(fil==2){
      result<-paste(result,fila1[fil])
    }else if(fil!=1) {
      result<-paste(result,"+(",fila1[fil])
      for(numF in 3:fil){
        result<-paste(result,"*(x-",columna1[numF-2],")")
      }
      result<-paste(result,")")
    }
    
    
  }
  return (result)
}


intNewton_table <- function(x, y) {
  
  c <- length(x)
  f<-length(x)

  matrizNan<-pracma::zeros(c,c+1)
  for (col in 0:(c+1)) {
    
    for (fil in 1:c ){
      matrizNan[fil,col]<-NaN
    }
  }
  
  
  z <- pracma::zeros(c, c+1)
  

  z[, 1] <- t(x)
  z[, 2] <- t(y)

  unicos<-unique(z[, 1])
  matrizNan[, 1] <- t(x)
  matrizNan[, 2] <- t(y)
  #Se comprueba que todos los datos introducidos no son 0 
  bool<-(x==pracma::zeros(c,1))
  cont<-0
  
      for(j in 1:c){
        if(bool[j,1]==TRUE){
         cont<-cont+1
        }
      }
  
  if(cont==(c)){
    return( pracma::zeros(c,1))
  }
 
  if( length(unicos)!=c){
    return( pracma::zeros(c,1))
  }
 
  
  for (col in 3:(c+1)) {
    
    for (fil in 1:(f-1)) {
      if((matrizNan[fil + 1, col - 1] - matrizNan[fil, col - 1])/(matrizNan[fil + col - 2] - matrizNan[fil])<0.000001 
         &&(matrizNan[fil + 1, col - 1] - matrizNan[fil, col - 1])/(matrizNan[fil + col - 2] - matrizNan[fil])>0
         ){
           
         matrizNan[fil, col] <-0                 
      
      }else if ((matrizNan[fil + 1, col - 1] - matrizNan[fil, col - 1])/(matrizNan[fil + col - 2] - matrizNan[fil])>-0.000001 
           &&(matrizNan[fil + 1, col - 1] - matrizNan[fil, col - 1])/(matrizNan[fil + col - 2] - matrizNan[fil])<0
        ){
          
          matrizNan[fil, col] <-0 
        }else{
      
          matrizNan[fil, col] <- round((matrizNan[fil + 1, col - 1] - matrizNan[fil,
                                              col - 1])/(matrizNan[fil + col - 2] - matrizNan[fil]),3)
      }
        if( matrizNan[fil, col]==Inf ||  matrizNan[fil, col]==-Inf|| is.nan(matrizNan[fil, col] )){
            matrizNan[fil, col]=NaN
           
        }
    }
   
    f<-f-1
    
  }
  return(matrizNan)
}


newtonII <- function(f, df, x0, n = 10, TOL = 1.e-5) {
  res <- data.frame(k=NULL,
                    x0 = NULL,
                    fx = NULL,
                    df = NULL,
                    x1=NULL)
  x <- x0
  
  for (k in seq(n)) {
    Fx <- f(x)
    
    DFx <- df(x)
    
    if (DFx == 0) {
      stop("Zero derivative.", call. = FALSE)
      
    }
    x1 <- x - (Fx / DFx)
    res <- rbind(res,
                 data.frame(
                   k=k-1,
                   x_k = x ,
                   `f(x)` = Fx,
                   `f'(x)`=DFx,
                   `x_k+1` = x1
                 ))
    x <- x - (Fx / DFx)
    
    
    if (is.nan(x) || is.infinite(x)) {
      stop("NaN", call. = FALSE)
    }
    
  }
  
  if (abs(f(x)) >= TOL) {
    warning("Tolerance not achieved.", call. = FALSE)
    
  }
  
  return(res)
}

#Ecuaciones no lineales
newton <- function(f, df, x0, n = 10, TOL = 1.e-5) {
  res <- data.frame(n = NULL,
                    x_n = NULL,
                    `f(x_n)` = NULL)
  x <- x0
  
  for (k in seq(n)) {
    Fx <- f(x)
    
    DFx <- df(x)
    
    if (DFx == 0) {
      stop("Zero derivative.", call. = FALSE)
      
    }
    x1 <- x - (Fx / DFx)
    res <- rbind(res,
                 data.frame(
                   n = k ,
                   x_n = x,
                   pto_corte=x1,
                   `f(x_n)` = Fx
                 ))
    x <- x - (Fx / DFx)
    
    
    if (is.nan(x) || is.infinite(x)) {
      stop("NaN", call. = FALSE)
    }
    
  }
  
  if (abs(f(x)) >= TOL) {
    warning("Tolerance not achieved.", call. = FALSE)
    
  }
  
  return(res)
}

bipart <- function(F, a, b, n) {
  res <- data.frame(a=NULL,b=NULL)
  FA <- feval(F, a)
  FB <- feval(F, b)
  
  x <- (a + b)/2
  
  e <- 0
  
  if (FA * FB > 0) {
    
    stop("There is no solution in the specified interval.",
         call. = FALSE)
  }
  
  if (FA == 0) {
    x <- a  # Si f(a) = 0, devolvemos x = a
    return(x)
  }
  if (FB == 0) {
    x <- b  # Si f(b) = 0, devolvemos x = b
    return(x)
  }
  
  if (n > 0) {
    
    
    for (k in 1:n) {
      
      FX <- feval(F, x)
      
      if (FX == 0) {
        return(x)
      }
      
      if (FA * FX < 0) {
        b <- x
      } else {
        a <- x
        FA <- feval(F, a)
      }  # (if)
      
      res <- rbind(res,
                   data.frame(a = a,
                              b = b))
      
      x <- (a + b)/2
      
    }  # (for k)
    
  }
  
  
  # tol <- (b - a)/2
  
  return(res)
  
}  # (function)



bipart_table <- function(F, a, b, n,k) {
  
  res <- data.frame(n = NULL, x_n = NULL,
                    `f(x_n)` = NULL,
                    I = NULL,
                    CotaError = NULL)
  
  FA <- feval(F, a)
  FB <- feval(F, b)
  
  x <- (a + b)/2
  
  e <- 0
  
  if (FA * FB > 0) {
    
    stop("There is no solution in the specified interval.",
         call. = FALSE)
  }
  
  if (FA == 0) {
    x <- a  # Si f(a) = 0, devolvemos x = a
    return(x)
  }
  if (FB == 0) {
    x <- b  # Si f(b) = 0, devolvemos x = b
    return(x)
  }
  
  if (n > 0)
    for (k in 1:n) {
      
      FX <- feval(F, x)
      cota <- 0.5 * (b - a)
      
      if (FX == 0) {
        return(x)
      }
      
      if (FA * FX < 0) {
        b <- x
      } else {
        a <- x
        FA <- feval(F, a)
      }  # (if)
      
      res <- rbind(res,
                   data.frame(n = k - 1,
                              x_n = x,
                              `f(x_n)` = FX,
                              Next_I = paste0("[", a, ", ", b, "]"),
                              CotaError = cota))
      
      x <- (a + b)/2
      
    }  # (for k)
  
  # tol <- (b - a)/2
  
  return(res)
  
}  # (function)


ui <- fluidPage(
  withMathJax(),
  navbarPage(
    theme = shinytheme("united"),
    "",
    ##MENU PRINCIPAL
    tabPanel("MÉTODOS NUMÉRICOS",
             mainPanel(
               fluidRow(column(
                 7,
                 br(),
                 h4(
                   "Esta web tiene como propósito ayudar a estudiantes de forma interactiva a resolver
        y entender diferentes ejercicios relacionados con contenidos de la asignatura MÉTODOS NUMÉRICOS. "
                 ),
                 br(),
                 h4(
                   "Puedes consultar información teórica y problemas resueltos sobre el temario que se encuentra en el menú superior.
                     Estos ejercicios se resuelven paso a paso, tanto los cálculos como las gráficas que se generan a partir de estos."
                 ),
                 h4(
                   "Encontrarás tanto ejemplos ya creados como otros donde seras tú el que deberá introducir ciertos datos del ejercicio para resolverlo."
                 ),
                 br(),
               ),
               br(),
               br()),
               
               fluidRow(column(
                 8,
                 br(),
                 h5("Trabajo fin de grado realizado por Delia Gómez Lobato."),
                 h5("Grado en Ingeniería Informática."),
                 br(),
                 column(
                   8,
                   img(
                     src = "rstudio.png",
                     height = 35,
                     width = 100,
                     align = "center"
                   ),
                   img(
                     src = "shiny.png",
                     height = 75,
                     width = 75,
                     align = "center"
                   )
                   
                 )
               ))
             )),
    
    
    
    
    ##APARTADO ECUACIONES LINEALES
    tabPanel("Ecuaciones no lineales",
             tabsetPanel(
               tabPanel("Teoría", mainPanel(fluidRow
                                            (
                                              column(
                                                6,
                                                br(),
                                                br(),
                                                h4("Método de bipartición: ", style =
                                                     "color: #f2540c"),
                                                br(),
                                                p("Dada una ecuación \\(f(x)=0\\) con \\(f(x)\\) continua."),
                                                
                                                column(
                                                  12,
                                                  offset = 0.5,
                                                  p("1) Buscar un intervalo \\([a , b ]\\) tal que \\(f(a)f(b)<0\\)"),
                                                  br(),
                                                  p("2) Repetir hasta criterio de parada:"),
                                                  column(
                                                    12,
                                                    offset = 0.5,
                                                    p("Tomar  \\(x=\\frac{a+b}{2}\\)."),
                                                    p(
                                                      "Evaluar \\(f(x):\\begin{cases}
                                               f(x)=0  & \\text{Parar,} \\\\
                                               f(a)f(x)<0 & \\text{Nuevo intervalo [a,x]}\\\\
                                               f(x)f(b)<0 & \\text{Nuevo intervalo [x,b]}
                                                       \\end{cases}\\! \\) "
                                                    ),
                                                    br(),
                                                  ),
                                                  p("3) Devolver el valor \\(x=\\frac{a+b}{2}\\) como solución. ")
                                                ),
                                              ),
                                              
                                              column(
                                                6,
                                                offset = 0.5,
                                                br(),
                                                br(),
                                                h4("Método de Newton: ", style = "color: #f2540c"),
                                                br(),
                                                p("Partir de un punto inicial \\(x_0\\) e iterar mediante: "),
                                                column(
                                                  12,
                                                  offset = 2,
                                                  p("\\(x_k+1= x_k -\\frac{f(x_k)}{f'(x_k)} \\)", style =
                                                      "font-size:large"),
                                                ),
                                              ),
                                            ))),
               
               
               ##BIPART
               tabPanel("Ejemplos Bipartición",
                        mainPanel(
                          navlistPanel(
                            tabPanel("Ecuación exponencial",  mainPanel(
                              fluidRow(
                                br(),
                                h4("Ecuación exponencial con intervalo [0,1]: $$f(x)=e^{(3x)}-4$$ "),
                                br(),
                              ),
                              fluidRow(
                                splitLayout(
                                  cellWidths = 500,
                                  fluidPage(verticalLayout(
                                    tableOutput("tablaBipartExp"), htmlOutput("textoBipartExp")
                                  )),
                                  
                                  plotOutput("graficaBipartExp", height = 600, width = 500)
                                ),
                                br(),
                                column(
                                  width = 12,
                                  offset = 8,
                                  align = "right",
                                  actionButton("botonBipartExp", ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                                ),
                                
                              ),
                              
                            )),
                            
                            tabPanel("Ecuación trigonométrica",  mainPanel(
                              fluidRow(
                                br(),
                                h4(
                                  "Ecuación trigonométrica con intervalo [-1,1]: $$f(x)=sin(x)+cos(x^2)$$ "
                                ),
                                br(),
                              ),
                              fluidRow(
                                splitLayout(
                                  cellWidths = 500,
                                  fluidPage(verticalLayout(
                                    tableOutput("tablaBipartTrig"),
                                    htmlOutput("textoBipartTrig")
                                  )),
                                  
                                  plotOutput("graficaBipartTrig", height = 600, width = 500)
                                ),
                                br(),
                                column(
                                  width = 12,
                                  offset = 8,
                                  align = "right",
                                  actionButton("botonBipartTrig", ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                                ),
                                
                              ),
                              
                            )),
                            
                            tabPanel(
                              "Punto de corte entre funciones",
                              mainPanel(
                                fluidRow(splitLayout(
                                  cellWidths = 500,
                                  fluidPage(verticalLayout(
                                    br(),
                                    h4(
                                      "Punto de corte entre las funciones $$f(x)=x^2-4x+5$$ $$g(x)=x+1$$ "
                                    ),
                                    h4("Ecuación con intervalo [2,7]: $$F(x)=x^2-5x+4$$ "),
                                    br(),
                                  )),
                                  fluidPage(verticalLayout(
                                    br(), br(), br(),
                                    h4("
               $$f(x)=g(x)$$ $$F(x)=f(x)-g(x)=0$$ ")
                                  ))
                                )),
                                fluidRow(
                                  splitLayout(
                                    cellWidths = 500,
                                    fluidPage(verticalLayout(
                                      tableOutput("tablaBipartPC"), htmlOutput("textoBipartPC")
                                    )),
                                    plotOutput("graficaBipartPC", height = 600, width = 500)
                                  ),
                                  br(),
                                  column(
                                    width = 12,
                                    offset = 8,
                                    align = "right",
                                    actionButton("botonBipartPC", ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                                  ),
                                  
                                ),
                                
                              )
                            ),
                            tabPanel("Optimización", mainPanel(
                              fluidRow(
                                br(),
                                h4(
                                  "Ejercicio de optimización con intervalo [0,2]: $$f(x)= (x^3-x-1)$$ $$f'(x)= (3x^2-1)$$"
                                ),
                                br(),
                              ),
                              fluidRow(
                                splitLayout(
                                  cellWidths = 500,
                                  fluidPage(verticalLayout(
                                    tableOutput("tablaBipartOpt"), htmlOutput("textoBipartOpt")
                                  )),
                                  plotOutput("graficaBipartOpt", height = 600, width = 500)
                                ),
                                br(),
                                column(
                                  width = 12,
                                  offset = 8,
                                  align = "right",
                                  actionButton("botonBipartOpt", ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                                ),
                                
                              ),
                              
                            )),
                            
                          )
                        )),
               
               ###NEWTON
               tabPanel(
                 "Ejemplos Newton",
                 uiOutput("Newton"),
                 mainPanel(navlistPanel(
                   tabPanel("Ecuación exponencial",  mainPanel(
                     fluidRow(br(),
                              splitLayout(
                                cellWidths = 500,
                                h4(
                                  " $$f(x)=e^{(x)}+x^2-4$$ $$f'(x)=e^{(x)}+2x$$ con \\(x_0=3.5\\)"
                                ),
                                fluidPage(verticalLayout(
                                  br(),
                                  br(),
                                  p("Recordamos:"),
                                  p("\\(x_{k+1}= x_k -\\frac{f(x_k)}{f'(x_k)} \\)", style = "font-size:large")
                                )),
                                
                              ), br()),
                     fluidRow(
                       splitLayout(
                         cellWidths = 500,
                         fluidPage(verticalLayout(
                           tableOutput("tablaNewtonExp"), htmlOutput("textoNewtonExp")
                         )),
                         
                         plotOutput("graficaNewtonExp", height = 600, width = 500),
                         
                       ),
                       br(),
                       column(
                         width = 12,
                         offset = 8,
                         align = "right",
                         actionButton("botonNewtonExp", ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                       ),
                       
                     ),
                     
                   )),
                   
                   tabPanel("Ecuación trigonométrica", mainPanel(
                     fluidRow(br(),
                              splitLayout(
                                cellWidths = 500,
                                h4("$$f(x)=x^2 -cos(x)$$$$f'(x)=2x+sen(x)$$ con \\(x_0=2\\) "),
                                fluidPage(verticalLayout(
                                  br(),
                                  br(),
                                  p("Recordamos:"),
                                  p("\\(x_{k+1}= x_k -\\frac{f(x_k)}{f'(x_k)} \\)", style = "font-size:large")
                                )),
                                
                              ), br()),
                     fluidRow(
                       splitLayout(
                         cellWidths = 500,
                         fluidPage(verticalLayout(
                           tableOutput("tablaNewtonTrig"),
                           htmlOutput("textoNewtonTrig")
                         )),
                         
                         plotOutput("graficaNewtonTrig", height = 600, width = 500),
                         
                       ),
                       br(),
                       column(
                         width = 12,
                         offset = 8,
                         align = "right",
                         actionButton("botonNewtonTrig", ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                       ),
                       
                     ),
                     
                   )),
                   
                   tabPanel(
                     "Punto de corte entre funciones",
                     mainPanel(
                       fluidRow(br(),
                                splitLayout(
                                  cellWidths = 500,
                                  h4(
                                    " $$f(x)=x^4+3x^2-10$$ $$g(x)=sen(x)-x^3$$ $$F(x)=x^4-x^3+3x^2-sen(x)-10$$ $$F'(x)=4x^3-3x^2+6x-cos(x)$$Con \\(x_0=-0.5\\) "
                                  ),
                                  fluidPage(verticalLayout(
                                    br(),
                                    br(),
                                    h4(" $$f(x)=g(x)$$ $$F(x)=f(x)-g(x)=0$$ "),
                                    p("Recordamos:"),
                                    p("\\(x_{k+1}= x_k -\\frac{f(x_k)}{f'(x_k)} \\)", style = "font-size:large")
                                  )),
                                  
                                ), br()),
                       fluidRow(
                         splitLayout(
                           cellWidths = 500,
                           fluidPage(verticalLayout(
                             tableOutput("tablaNewtonPC"), htmlOutput("textoNewtonPC")
                           )),
                           plotOutput("graficaNewtonPC", height = 600, width = 500)
                         ),
                         br(),
                         column(
                           width = 12,
                           offset = 8,
                           align = "right",
                           actionButton("botonNewtonPC", ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                         ),
                         
                       ),
                       
                     )
                   ),
                   
                   tabPanel("Optimización",  mainPanel(
                     fluidRow(br(),
                              splitLayout(
                                cellWidths = 500,
                                h4(
                                  "Ejercicio optimización $$f(x)=x^3+4x^2-10$$ $$f'(x)=3x^2+8x$$ $$f''(x)=6x+8$$ con \\(x_0=0.75\\) "
                                ),
                                fluidPage(verticalLayout(
                                  br(),
                                  br(),
                                  p("En este caso:"),
                                  p("\\(x_{k+1}= x_k -\\frac{f'(x_k)}{f''(x_k)} \\)", style = "font-size:large")
                                )),
                                
                              ), br()),
                     fluidRow(
                       splitLayout(
                         cellWidths = 500,
                         fluidPage(verticalLayout(
                           tableOutput("tablaNewtonOpt"), htmlOutput("textoNewtonOpt")
                         )),
                         
                         plotOutput("graficaNewtonOpt", height = 600, width = 500),
                         
                       ),
                       br(),
                       column(
                         width = 12,
                         offset = 8,
                         align = "right",
                         actionButton("botonNewtonOpt", ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                       ),
                       
                     ),
                     
                   )),
                   
                 )),
               )
             )), 
    
    ##APARTADO INTERPOLACIÓN
    navbarMenu(
      "Interpolación",
      tabPanel("Polinómica",
               mainPanel(fluidRow
                         (
                           column(
                             12,
                             br(),
                             h4("Problema clásico: "),
                             br(),
                             p(
                               "Consiste en hallar una función polinómica del menor grado posible que pase por los puntos dados."
                             ),
                             
                             column(
                               12,
                               offset = 0.5,
                               p(
                                 "Si consideramos la base {\\(1,x,x^2,...,x^n\\)}  el polinomio \\(P_n (x)\\) que buscamos es
             $$P_n (x)= a_0+a_1 x+a_2 x^2+...+a_n x^n$$"
                               ),
                               
                               p(
                                 "Como queremos que \\(P_n(x_i) = f(x_i)\\) con \\(i =0,1,...,n\\), para encontrar los coeficientes
                \\(a_0,a_1,a_2,...,a_n\\) hemos de resolver el sistema de ecuaciones lineales :"
                               ),
                               column(
                                 12,
                                 offset = 0.5,
                                 p(
                                   "$$\\left. \\begin{align}
                  \\  a_0+a_1 x_0+a_2 x_0^2+...+a_n x_0^n &= f(x_0) \\\\
                 \\  a_0+a_1 x_1+a_2 x_1^2+...+a_n x_1^n &= f(x_1) \\\\
                 \\ . \\\\
                 \\ . \\\\
                 \\ . \\\\
                 \\ a_0+a_1 x_n+a_2 x_n^2+...+a_n x_n^n &= f(x_n)
                 \\end{align}\\right\\} \\Rightarrow
                 \\begin{pmatrix}
                    1 & x_0 & x_0^2 & . &.&.& x_0^n\\\\
                    1 & x_1 & x_1^2 & . &.&.& x_1^n \\\\
                    . & .   & .     & . & & & . \\\\
                    . & .   & .     &  &.& & .\\\\
                    . & .   & .     &  & & .& .\\\\
                    1 & x_n & x_n^2 & . &.&.& x_n^n
                    \\end{pmatrix}
               \\begin{pmatrix}
                   a_0\\\\
                   a_1\\\\
                    .  \\\\
                    .   \\\\
                    .    \\\\
                    a_n
                    \\end{pmatrix} =
                    \\begin{pmatrix}
                   f(x_0)\\\\
                   f(x_1)\\\\
                    .  \\\\
                    .   \\\\
                    .    \\\\
                    f(x_n)
                    \\end{pmatrix} $$"
                                 )
                               ),
                               br(),
                               br(),
                               p(
                                 "La matriz de coeficientes anterior es una matriz de Vandermonde, cuyo determinante viene dado por \\( \\prod\\nolimits_{j<k}(x_j-x_k)\\) .
       Por tanto, si todos los nodos son distintos, dicho determinante es distinto de
        cero, y el sistema anterior tiene una única solución. "
                               ),
                               br(),
                               
                               p(
                                 "El error cometido al aproximar\\(f(x)\\) por \\(P_n(x)\\) se puede expresar como  "
                               ),
                               
                               p(
                                 "$$ e_n(x) = \\frac{f^{n+1}(\\xi)}{(n+1)!}\\ (x-x_0)(x-x_1)...(x-x_n) $$ "
                               ),
                               p("con \\(\\xi\\) entre el mayor y el menor de los nodos. ")
                             ),
                           ),
                         )),),
      
      tabPanel(
        "De Newton",
        tabsetPanel(
          tabPanel("Teoría", 
                   mainPanel(fluidRow
                             (
                               column(
                                 12,
                                 br(),
                                 h4("Interpolación de Newton en diferencias divididas"),
                                 br(),
                                 p("Se llama diferencias divididas a los coeficientes ."),
                                 p("$$f [ x_j,x_j+1,...,x_i-1,x_i]$$"),
                                 p("Se pueden calcular usando el esquema:"),
                                 column(
                                   12,
                                   offset = 0.5,
                                   p(
                                     "$$\\begin{array}{|c|c|}
      i & x_i & f(x_i) & f[x_i] & f[x_i-1,x_i] & f[x_i-2,x_i-1,x_i] &f[x_i-3,x_i-2,x_i-1,x_i] \\\\
      \\hline
      0 & x_0  & f(x_0)  & f[x_0] \\\\
                                     &&&& f[x_0,x_1]=\\frac{f(x_1)-f(x_0)}{x_1-x_0}\\\\
      1 & x_1 & f(x_1)  & f[x_1]&& f[x_0,x_1,x_2]=\\frac{f[x_1,x_2]-f[x_0,x_1]}{x_2-x_0}\\\\
                                     &&&& f[x_1,x_2]=\\frac{f(x_2)-f(x_1)}{x_2-x_1}&& f[x_0,x_1,x_2,x_3]\\\\
      2 & x_2& f(x_2)  & f[x_2]&& f[x_1,x_2,x_3]=\\frac{f[x_2,x_3]-f[x_1,x_2]}{x_3-x_1} \\\\
                                     &&&& f[x_2,x_3]=\\frac{f(x_3)-f(x_2)}{x_3-x_2}\\\\
      3 & x_3& f(x_3)  & f[x_3]
      \\end{array}$$"
                                   )
                                 ),
                                 p("Aunque en algunos casos también se expresa de la siguiente manera:"),
                                 column(
                                   12,
                                   offset = 0.5,
                                   p(
                                     "$$\\begin{array}{|c|c|}
      i & x_i & f(x_i) & f[x_i] & f[x_i-1,x_i] & f[x_i-2,x_i-1,x_i] &f[x_i-3,x_i-2,x_i-1,x_i] \\\\
      \\hline
      0 & x_0  & f(x_0)  & f[x_0] & f[x_0,x_1]& f[x_0,x_1,x_2]& f[x_0,x_1,x_2,x_3]\\\\
                                     &&&\\\\
      1 & x_1 & f(x_1)  & f[x_1]& f[x_1,x_2]& f[x_1,x_2,x_3]\\\\
                                     &&&\\\\
      2 & x_2& f(x_2)  & f[x_2]&f[x_2,x_3]\\\\
                                     &&&& \\\\
      3 & x_3& f(x_3)  & f[x_3]
      \\end{array}$$"
                                   )
                                 ),
                                 p(
                                   "En la primera fila de la tabla se observan los coeficientes del polinomio de interpolación."
                                 ),
                                 br(),
                                 strong("Fórmula de Newton", style = "color: #f2540c"),
                                 p(
                                   "El polinomio de interpolación se pude escribir como: $$ P_n(x) = f[x_0]+ f[x_0,x_1](x-x_0)+f[x_0,x_1,x_2](x-x_0)(x-x_1)+...+f[x_0,x_1,x_2,...,x_n](x-x_0)(x-x_1)...(x-x_n).$$"
                                 ),
                                 strong("Error de interpolación", style = "color: #f2540c"),
                                 p(
                                   "Si \\(f\\) es una función continua en \\([a , b ], K \\) veces derivables en \\((a , b )\\) y la derivada k-ésima es continua, entonces:
$$ E(x) = \\frac{f^{n+1}(\\xi)}{(n+1)!}\\prod\\limits_{i=0}^n (x-x_i) $$ "
                                   
                                 ),
                               )
                             ))),
          tabPanel("Ejemplos", mainPanel(
            fluidRow(
              br(),
              h4("Interpolación de Newton en diferencias divididas", style =
                   "color: #f2540c"),
              br(),
            ),
            fluidRow(
              splitLayout(
                cellWidths = 800,
                DT::dataTableOutput("InterpolacionNewton"),
                plotOutput("graficaInterpolacionNewton")
              ),
             
              htmlOutput("textoInterpolacionNewton"),
              br(),
              fluidRow(
                splitLayout(
                  cellWidths = 800,
                  fluidPage(
                    verticalLayout(
                      p("Introduce los puntos y pulsa siguiente paso cuando hayas terminado : "),
                      numericInput("x", "X:", 0, min = -50, max = 50),
                      numericInput("y", "Y:", 0, min = -50, max = 50),
                      splitLayout(
                        cellWidths = 150,
                        actionButton("botonNewtonPuntos", "Introducir punto", style = "color: #fff; background-color: #f2540c; border-color: #f2540c"),
                        actionButton("botonInterpolacionNewton", "Sig. Paso", style = "color: #fff; background-color: #f2540c; border-color: #f2540c"),
                        actionButton("botonNewtonReset", "Reset", style = "color: #fff; background-color: #f2540c; border-color: #f2540c"),))),
                   
                ),
                
                
                
              ),
              
            )),
          ))
      ),
      tabPanel(
        "Inversa",
      
        tabsetPanel(
          tabPanel("Teoría",
                   mainPanel(fluidRow
                             (
                               column(
                                 12,
                                 br(),
                                 h4("Interpolación inversa"),
                                 br(),
                                 
                                 p(
                                   "Dados unos valores de la función \\(\\left\\{(x_i,y_i)\\right\\}\\) y un valor \\(y\\), el objetivo ahora es encontrar el valor de \\(x\\), tal que \\(f(x)=y\\)  ."
                                 ),
                                 column(
                                   12,
                                   offset = 0.5,
                                   br(),
                                   strong("Método 1:", style = "color: #f2540c"),
                                   p(
                                     "Si los valores \\(y_i\\) son diferentes entre sí, podemos cambiar los papeles de \\(x_i\\) e \\(y_i\\). Lo que hacemos es interpolar a los datos \\((x_i,y_i)\\). "
                                   ),
                                   br(),
                                   strong("Método 2:", style = "color: #f2540c"),
                                   p(
                                     "Interpolamos una función a los datos \\(\\left\\{(x_i,y_i)\\right\\}\\), de manera que se obtiene una estimación de \\(f(x)\\).Después resolvemos la ecuación \\(f(x)=y\\)  . "
                                   )
                                 ),
                               )
                             ))),
          tabPanel("Ejemplos", mainPanel(
            fluidRow(
              br(),
              h4("Interpolación inversa", style =
                   "color: #f2540c"),
              br(),
            ),
            fluidRow(
              splitLayout(
                cellWidths = 800,
                DT::dataTableOutput("InterpolacionInversa"),
                plotOutput("graficaInterpolacionInversa")
              ),
              htmlOutput("textoInterpolacionInversa"),
              br(),
              fluidRow(
                splitLayout(
                  cellWidths = 500,
                  fluidPage(
                    verticalLayout(
                      p("Introduce los puntos y pulsa siguiente paso cuando hayas terminado : "),
                      
                      numericInput("x1", "X:", 0, min = -50, max = 50),
                      numericInput("y1", "Y:", 0, min = -50, max = 50),
                      
                      splitLayout(
                        cellWidths = 150,
                        actionButton("botonInversaPuntos", "Introducir punto", style = "color: #fff; background-color: #f2540c; border-color: #f2540c"),
                        actionButton("botonInterpolacionInversa", "Sig. Paso", style = "color: #fff; background-color: #f2540c; border-color: #f2540c"),
                        actionButton("botonInversaReset", "Reset", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")),
                      br(),)),
                  fluidPage(
                    verticalLayout(
                      htmlOutput("textoInterpolacionValor"),
                      br(),
                      uiOutput("y2"),
                      uiOutput("botonIntInvV"))))),
                
                ),
                
                
                
              ),
              
            )),
          ),
      tabPanel(
        "Osculatoria",
        tabsetPanel(
          tabPanel("Teoría", 
                   mainPanel(fluidRow
                             (
                               column(
                                 12,
                                 br(),
                                 h4("Interpolación osculatoria"),
                                 br(),
                                 
                                 p(
                                   "Este método se utiliza si además de exigir que el polinomio interpolador coincida con la función en los nodos, se impone que ciertas derivadas en dichos nodos también coincidan con las
correspondientes deirvadas de la función."
                                 ),
                                 column(
                                   12,
                                   offset = 0.5,
                                   
                                   p(
                                     "Utilizaremos el esquema en diferencias divididas visto en la interpolación de Newton, pero repitiendo \\(k_i\\) veces cada nodo \\(x_i\\) y utilizando que:
$$ f[x_i,x_i]=f'(x_i), f[x_i,x_i,x_i]=\\frac{f''(x_i)}{2},...,f[x_i,x_i,...^{n+1},x_i]=\\frac{f^{n}(x_i)}{n!} $$"
                                   ),
                                   br(),
                                   
                                 ),
                                 strong("Interpolación de Hermite", style = "color: #f2540c"),
                                 p(
                                   "Hablamos de este caso particular cuando se usa como información el valor de la función y de su derivada en cada uno de los nodos.
    Se trata del polinomio de interpolación osculatoria con \\(k_i=1, i = 0,...,n\\), y de grado por tanto \\(2n+1\\). "
                                 ),
                                 p(
                                   "La fórmula de error en este caso es: $$ e_{2n+1}(x) = \\frac{f^{2n+2}(\\xi)}{(2n+2)!}\\ (x-x_0)^2(x-x_1)^2...(x-x_n)^2$$"
                                 )
                               )
                             ))),
          tabPanel("Ejemplos", mainPanel(
            fluidRow(
              br(),
              h4("Interpolación de Hermite", style =
                   "color: #f2540c"),
              br(),
            ),
            fluidRow(
              splitLayout(
                cellWidths = 800,
                DT::dataTableOutput("InterpolacionHermite"),
                plotOutput("graficaInterpolacionHermite")
              ),
              htmlOutput("textoInterpolacionHermite"),
              br(),
              fluidRow(
                splitLayout(
                  cellWidths = 500,
                  fluidPage(
                    verticalLayout(
                      p("Introduce los puntos y pulsa siguiente paso cuando hayas terminado : "),
                      numericInput("xh", "x:", 0, min = -50, max = 50),
                      numericInput("yh", "f(x):", 0, min = -50, max = 50),
                      numericInput("yhd", "f(x)':", 0, min = -50, max = 50),
                      splitLayout(
                        cellWidths = 150,
                        actionButton("botonHermitePuntos", "Introducir datos", style = "color: #fff; background-color: #f2540c; border-color: #f2540c"),
                        actionButton("botonInterpolacionHermite", "Sig. Paso", style = "color: #fff; background-color: #f2540c; border-color: #f2540c"),
                        actionButton("botonHermiteReset", "Reset", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")),
                      br(),)),
              
            ),
            
            
            
          ),
          
          )),
      ),
        
        ))),
    
    #APARTADO DERIVACION
    tabPanel("Derivación",
             
             tabsetPanel(
               tabPanel("Teoría",
                        mainPanel(fluidRow(
                          column(
                            12,
                            br(),
                            br(),
                            h4("Fórmulas de derivación numérica "),
                            br(),
                            strong(" Fórmula de dos puntos centrada", style = "color: #f2540c"),
                            br(),
                            p(
                              "$$\\begin{array}{c|cc}
xi & ... &  x_0 - h &  x_0 & x_0 + h & ... \\\\
\\hline
f(x_i) & ...  & f(x_0 - h)  & f(x_0) & f(x_0 + h) & ... \\\\
\\end{array}$$",
                              style = "font-size:small"
                            ),
                            
                            p(
                              "Teniendo en cuenta que  \\((\\frac{f'''(\\xi)}{3!}-h^2) \\) representa el error y que (\\( x_0-h<\\xi< x_0+h\\)):"
                            ),
                            br(),
                            column(
                              11,
                              offset = 4,
                              p(
                                "\\(f'(x_0) = \\frac{f(x_0+h)-f(x_0-h)}{2h}-(\\frac{f'''(\\xi)}{3!}-h^2) \\)",
                                style = "font-size:medium"
                              ),
                              br(),
                            ),
                            p("Es la más fiable de las tres ya que da lugar a un error menor."),
                            br(),
                            strong(" Fórmula de dos puntos descentrada a la derecha", style =
                                     "color: #f2540c"),
                            br(),
                            p(
                              "$$\\begin{array}{c|cc}
xi & x_0 & x_0 + h & ... \\\\
\\hline
f(x_i) & f(x_0) & f(x_0 + h) & ... \\\\
\\end{array}$$",
                              style = "font-size:small"
                            ),
                            p(
                              "Teniendo en cuenta que  \\((\\frac{f'''(\\xi)}{2!}-h) \\) representa el error y que (\\( x_0<\\xi< x_0+h\\))"
                            ),
                            br(),
                            column(
                              11,
                              offset = 4,
                              p(
                                "\\(f'(x_0) = \\frac{f(x_0+h)-f(x_0)}{h}-(\\frac{f'''(\\xi)}{2!}-h) \\)",
                                style = "font-size:medium"
                              ),
                              br(),
                            ),
                            br(),
                            
                            strong("Fórmula de dos puntos descentrada a la izquierda", style =
                                     "color: #f2540c"),
                            br(),
                            p(
                              "$$\\begin{array}{c|cc}
xi & ... &  x_0 - h &  x_0  \\\\
\\hline
f(x_i) & ...  & f(x_0 - h)  & f(x_0) \\\\
\\end{array}$$",
                              style = "font-size:small"
                            ),
                            p(
                              "Teniendo en cuenta que  \\((\\frac{f'''(\\xi)}{2!}-h) \\)  representa el error y que (\\(x_0-h<\\xi< x_0\\))"
                            ),
                            br(),
                            column(
                              11,
                              offset = 4,
                              p(
                                "\\(f'(x_0) = \\frac{f(x_0)-f(x_0-h)}{h}-(\\frac{f'''(\\xi)}{2!}-h) \\)",
                                style = "font-size:medium"
                              ),
                              br(),
                            ),
                            br(),
                          )
                        ))),
               tabPanel(
                 "Ejemplos",
                 mainPanel(navlistPanel(
                   tabPanel("Ejemplo 1 ", mainPanel(
                     fluidRow(
                       splitLayout(
                         cellWidths = 325,
                         h4(""),
                         h4(
                           " $$f(x)=x^{(e-x^2)}e^{(-x^2)}(-2xlnx+\\frac{1}{x})$$  $$x=1.3$$ "
                         ),
                         h4("")
                       ),
                       br(),
                       splitLayout(
                         cellWidths = 500,
                         fluidPage(verticalLayout(
                           radioButtons(
                             "metodosDerivacion1",
                             "Selecciona una fórmula:",
                             c("Centrada", "Descentrada izquierda", "Descentrada derecha"),
                             selected = "Centrada",
                             inline = TRUE
                           ),
                           
                         )),
                         fluidPage(verticalLayout(
                           sliderInput(
                             "h",
                             "h:",
                             min = 0.01,
                             max = 0.2,
                             value = 0.05
                           )
                         )),
                         
                       ),
                       br()
                     ),
                     fluidRow(splitLayout(
                       cellWidths = 500,
                       fluidPage(verticalLayout(uiOutput("textoDerivacion1"),)),
                       plotOutput("graficaDerivacion1", height = 600, width = 500),
                       
                     ),),
                     
                   )),
                   tabPanel("Ejemplo 2", mainPanel(
                     fluidRow(
                       splitLayout(
                         cellWidths = 325,
                         h4(""),
                         h4(" $$f(x)=(x^2+3x-2)^4$$  $$x=1.3$$ "),
                         h4("")
                       ),
                       br(),
                       splitLayout(
                         cellWidths = 500,
                         fluidPage(verticalLayout(
                           radioButtons(
                             "metodosDerivacion2",
                             "Selecciona una fórmula:",
                             c("Centrada", "Descentrada izquierda", "Descentrada derecha"),
                             selected = "Centrada",
                             inline = TRUE
                           ),
                           
                         )),
                         fluidPage(verticalLayout(
                           sliderInput(
                             "h2",
                             "h:",
                             min = 0.01,
                             max = 0.2,
                             value = 0.05
                           )
                         )),
                         
                       ),
                       br()
                     ),
                     fluidRow(splitLayout(
                       cellWidths = 500,
                       fluidPage(verticalLayout(
                         uiOutput("textoDerivacion2"),
                         uiOutput("analiticoDerivacion"),
                       )),
                       plotOutput("graficaDerivacion2", height = 600, width = 500),
                       
                     ),),
                     
                   )),
                   
                 ))
               )
             )),
    
    #APARTADO INTEGRACIÓN
    navbarMenu(
      "Integración",
      tabPanel(
        "Fórmula cerrada de Newton-Côtes",
        tabsetPanel(
          tabPanel("Teoría", 
                   mainPanel(fluidRow(
                     column(
                       12,
                       br(),
                       br(),
                       h4("Fórmula cerrada de Newton-Côtes "),
                       br(),
                       p("La fórmula cerrada de Newton-Côtes de n+1 puntos utiliza los nodos: "),
                       p(
                         "\\(x_i=x_0+ih\\) , \\( i = 0,...,n \\) con \\(x_0=a\\) , \\( x_n =b\\) , y \\(h = \\frac{b-a}{n} \\)"),
                       p("Su grado de exactitud es n+1 si n par, y n si n impar."),
                       strong(" Regla del trapecio (n=1)", style =
                                "color: #f2540c"),
                       br(),br(),
                       p(
                         "Teniendo en cuenta que  \\((\\frac{h^3}{12}f''(\\xi)) \\) representa el error y que (\\( x_0<\\xi< x_1\\)):"
                       ),
                       
                       
                       br(),
                       column(
                         11,
                         offset = 4, 
                         p(
                           "\\(\\int_{x_0}^{x_1}f(x)dx =\\frac{h}{2}(f(x_0)+f(x_1)) - \\frac{h^3}{12}f''(\\xi)\\)",
                           style = "font-size:medium"
                         ),
                         
                         br(),
                       ),
                       
                       br(),
                       strong(" Regla de Simpson (n=2)", style =
                                "color: #f2540c"),
                       br(),br(),
                       p(
                         "Teniendo en cuenta que  \\((\\frac{h^5}{90}f^{iv}(\\xi)) \\) representa el error y que (\\( x_0<\\xi< x_2\\)):"
                       ),
                       br(),
                       column(
                         11,
                         offset = 4, 
                         p(
                           "\\(\\int_{x_0}^{x_2}f(x)dx =\\frac{h}{3}(f(x_0)+4f(x_1)+f(x_2)) - \\frac{h^5}{90}f^{iv}(\\xi)\\)",
                           style = "font-size:medium"
                         ),
                         
                         br(),
                       ),
                       
                       br(),
                       
                       strong(" Regla de Simpson tres octavos (n=3)", style =
                                "color: #f2540c"),
                       
                       br(),br(),
                       p(
                         "Teniendo en cuenta que  \\((\\frac{3h^5}{80}f^{iv}(\\xi)) \\) representa el error y que (\\( x_0<\\xi< x_3\\)):"
                       ),
                       br(),
                       column(
                         11,
                         offset = 4, 
                         p(
                           "\\(\\int_{x_0}^{x_3}f(x)dx =\\frac{3h}{8}(f(x_0)+3f(x_1)+3f(x_2)+f(x_3)) - \\frac{3h^5}{80}f^{iv}(\\xi)\\)",
                           style = "font-size:medium"
                         ),
                         
                         br(),
                       ),
                       
                       br(),
                       
                     )
                   ))),
          tabPanel(
            "Ejemplos",
            mainPanel(navlistPanel(
              
              tabPanel("Ejemplo 1 ",mainPanel(
                h4("Fórmula cerrada de Newton-Côtes", style =
                     "color: #f2540c"),
                fluidRow(
                  splitLayout(cellWidths = 325,h4(""),
                              h4("  $$\\int_{2}^{4}\\frac{x-2}{\\sqrt{(x-1)}}dx  $$ "), h4("")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasIntegracionNewtonC1", "Selecciona un método:",
                                    c("Regla del trapecio", "Regla de Simpson", "Regla de Simpson tres octavos"),
                                    selected = "Regla del trapecio",
                                    inline = TRUE
                                  ),
                                  
                                )),
                              
                              
                              
                  ),br()),
                fluidRow(
                  splitLayout(
                    cellWidths = 600,
                    fluidPage( 
                      verticalLayout(
                        
                          uiOutput("DatosIntegracionNewtonC1"),
                          uiOutput("TablaIntegracionNewtonC1"),
                          uiOutput("textIntegracionNewtonC1"),
                          uiOutput("resultIntegracionNewtonC1"),
                          
                          uiOutput("botonTrapecioNewtonC1"),
                          uiOutput("botonSimpsonNewtonC1"),
                          uiOutput("botonSimpson3NewtonC1")
                        )),
                    plotOutput("graficaIntegracionNewtonC1",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
              tabPanel("Ejemplo 2",mainPanel(
                h4("Fórmula cerrada de Newton-Côtes", style =
                     "color: #f2540c"),
                fluidRow(
                  splitLayout(cellWidths = 325,h4(""),
                              h4(" $$\\int_{0}^{4}x\\sqrt{(x^2+9)}dx $$ "), h4("")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasIntegracionNewtonC2", "Selecciona un método:",
                                    c("Regla del trapecio", "Regla de Simpson", "Regla de Simpson tres octavos"),
                                    selected = "Regla del trapecio",
                                    inline = TRUE
                                  ),
                                  
                                )),
                              
                              
                  ),br()),
                fluidRow(
                  splitLayout(
                    cellWidths = 600,
                    fluidPage( 
                      verticalLayout(
                        uiOutput("DatosIntegracionNewtonC2"),
                        uiOutput("TablaIntegracionNewtonC2"),
                        uiOutput("textIntegracionNewtonC2"),
                        uiOutput("resultIntegracionNewtonC2"),
                        uiOutput("AnaliticoIntegracionNewtonC2"),
                        
                        uiOutput("botonTrapecioNewtonC2"),
                        uiOutput("botonSimpsonNewtonC2"),
                        uiOutput("botonSimpson3NewtonC2")
                      )),
                    plotOutput("graficaIntegracionNewtonC2",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
              
            ))
          )
        )),
      
      
      
      tabPanel(
        "Fórmula abierta de Newton-Côtes",
        tabsetPanel(
          tabPanel("Teoría",
                   mainPanel(fluidRow(
                     column(
                       12,
                       br(),
                       br(),
                       h4("Fórmula abierta de Newton-Côtes "),
                       br(),
                       p("La fórmula abierta de Newton-Côtes de n+1 puntos utiliza los nodos: "),
                       p(
                         "\\(x_i=x_0+ih\\) con \\( i = 0,...,n \\) , \\(x_0=a+h\\) , \\( x_n =b-h(x_{-1}=a,x_{n+1}=b)\\) , y \\(h = \\frac{b-a}{n+2} \\)"),
                       p("Su grado de exactitud es n+1 si n par, y n si n impar."),
                       
                       strong(" Regla del punto medio (n=0)", style =
                                "color: #f2540c"),
                       br(),br(),
                       p(
                         "Teniendo en cuenta que  \\((\\frac{h^3}{3}f''(\\xi)) \\) representa el error y que (\\( x_{-1}<\\xi< x_1\\)):"
                       ),
                       
                       
                       br(),
                       column(
                         11,
                         offset = 4, 
                         p(
                           "\\(\\int_{x_{-1}}^{x_1}f(x)dx =2h(f(x_0)) + \\frac{h^3}{3}f''(\\xi)\\)",
                           style = "font-size:medium"
                         ),
                         
                         br(),
                       ),
                       
                       br(),
                       strong(" Fórmula abierta de dos puntos (n=1)", style =
                                "color: #f2540c"),
                       br(),br(),
                       p(
                         "Teniendo en cuenta que  \\((\\frac{3h^3}{4}f''(\\xi)) \\) representa el error y que (\\( x_{-1}<\\xi< x_2\\)):"
                       ),
                       br(),
                       column(
                         11,
                         offset = 4, 
                         p(
                           "\\(\\int_{x_{-1}}^{x_2}f(x)dx =\\frac{3h}{2}(f(x_0)+f(x_1)) + \\frac{3h^3}{4}f''(\\xi) \\)",
                           style = "font-size:medium"
                         ),
                         
                         br(),
                       ),
                       
                       br(),
                       
                       strong(" Fórmula abierta de tres puntos (n=2)", style =
                                "color: #f2540c"),
                       
                       br(),br(),
                       p(
                         "Teniendo en cuenta que  \\((\\frac{14h^5}{45}f^{iv}(\\xi) )\\) representa el error y que (\\( x_{-1}<\\xi< x_3\\)):"
                       ),
                       br(),
                       column(
                         11,
                         offset = 4, 
                         p(
                           "\\(\\int_{x_{-1}}^{x_3}f(x)dx =\\frac{4h}{3}(2f(x_0)-f(x_1)+2f(x_2)) + \\frac{14h^5}{45}f^{iv}(\\xi)\\)",
                           style = "font-size:medium"
                         ),
                         
                         br(),
                       ),
                       
                       br(),
                       
                     )
                   ))),
          tabPanel(
            "Ejemplos",
            mainPanel(navlistPanel(
              tabPanel("Ejemplo 1 ",mainPanel(
                h4("Fórmula abierta de Newton-Côtes", style =
                     "color: #f2540c"),
                fluidRow(
                  splitLayout(cellWidths = 325,h4(""),
                              h4(" $$\\int_{3}^{5}\\frac{x^3}{(x^2+1)(x^2+4)(x^2+9)}dx $$ "), h4("")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasIntegracionNewtonA1", "Selecciona un método:",
                                    c("Regla del punto medio", "Fórmula de 2 puntos", "Fórmula de 3 puntos"),
                                    selected = "Regla del punto medio",
                                    inline = TRUE
                                  ),
                                  
                                )),
                              
                              
                              
                  ),br()),
                fluidRow(
                  splitLayout(
                    cellWidths = 600,
                    fluidPage( 
                      verticalLayout(
                        uiOutput("DatosIntegracionNewtonA1"),
                        uiOutput("TablaIntegracionNewtonA1"),
                        uiOutput("textIntegracionNewtonA1"),
                        uiOutput("resultIntegracionNewtonA1"),
                        
                        uiOutput("botonPMNewtonA1"),
                        uiOutput("boton2PNewtonA1"),
                        uiOutput("boton3PNewtonA1")
                        )),
                    plotOutput("graficaIntegracionNewtonA1",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
              tabPanel("Ejemplo 2 ",mainPanel(
                h4("Fórmula abierta de Newton-Côtes", style =
                     "color: #f2540c"),
                fluidRow(
                  splitLayout(cellWidths = 325,h4(""),
                              h4("  $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx $$  "), h4("")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasIntegracionNewtonA2", "Selecciona un método:",
                                    c("Regla del punto medio", "Fórmula de 2 puntos", "Fórmula de 3 puntos"),
                                    selected = "Regla del punto medio",
                                    inline = TRUE
                                  ),
                                  
                                )),
                              
                              
                              
                  ),br()),
                
                fluidRow(
                  splitLayout(
                    cellWidths = 600,
                    
                    fluidPage( 
                      verticalLayout(
                        uiOutput("DatosIntegracionNewtonA2"),
                        uiOutput("TablaIntegracionNewtonA2"),
                        uiOutput("textIntegracionNewtonA2"),
                        uiOutput("resultIntegracionNewtonA2"),
                        uiOutput("AnaliticoIntegracionNewtonA2"),
                       
                          uiOutput("botonPMNewtonA2"),
                          uiOutput("boton2PNewtonA2"),
                          uiOutput("boton3PNewtonA2")
                        
                      )),
                    plotOutput("graficaIntegracionNewtonA2",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
            ))
          )
          
        )
      ),
      tabPanel(
        "Fórmulas de integración compuestas",
        uiOutput("ICompuesta"),tabsetPanel(
          tabPanel("Teoría",
                   mainPanel(fluidRow(
                     column(
                       12,
                       br(),
                       br(),
                       h4("Fórmulas de integración compuestas "),
                       br(),
                       p("Estas fórmulas resultan de dividir el intervalo de integración en subintervalos para aplicar una fórmula simple en cada uno de ellos, y así simplificar los cáculos. "),
                       
                       strong(" Regla del trapecio compuesta", style =
                                "color: #f2540c"),
                       br(),br(),
                       
                       p(
                         " Siendo \\(h=\\frac{b-a}{N}\\), teniendo en cuenta que  \\((\\frac{h^2(b-a)}{12}f''(\\xi)) \\) representa el error y que (\\( a<\\xi< b\\)):"
                       ),
                       
                       
                       br(),
                       column(
                         11,
                         offset = 4, 
                         p(
                           "\\(\\int_{a}^{b}f(x)dx =\\frac{h}{2}(f(a)+2\\sum\\limits_{j=1}^{N-1}f(x_j)+f(b)) - \\frac{h^2(b-a)}{12}f''(\\xi)\\)",
                           style = "font-size:medium"
                         ),
                         
                         br(),
                       ),
                       
                       br(),
                       strong(" Regla de Simpson compuesta", style =
                                "color: #f2540c"),
                       br(),br(),
                       p(
                         "Siendo \\(h=\\frac{b-a}{N}, N =2m\\), teniendo en cuenta que  \\((\\frac{h^4(b-a)}{180}f^{iv}(\\xi)) \\) representa el error y que (\\( a <\\xi< b\\)):"
                       ),
                       br(),
                       column(
                         11,
                         offset = 4, 
                         p(
                           "\\(\\int_{a}^{b}f(x)dx =\\frac{h}{3}(f(a)+2\\sum\\limits_{j=1}^{m-1}f(x_{2j})+4\\sum\\limits_{j=1}^{m}f(x_{2j-1})+f(b)) - \\frac{h^4(b-a)}{180}f^{iv}(\\xi)\\)",
                           style = "font-size:medium"
                         ),
                         
                         br(),
                       ),
                       
                       br(),
                       
                       strong("Regla del punto medio compuesta", style =
                                "color: #f2540c"),
                       
                       br(),br(),
                       p(
                         "Siendo \\(h=\\frac{b-a}{N+2}, N =2m\\), teniendo en cuenta que  \\((\\frac{h^2(b-a)}{6}f''(\\xi) )\\) representa el error y que (\\( a<\\xi<b\\)):"
                       ),
                       br(),
                       column(
                         11,
                         offset = 4, 
                         p(
                           "\\(\\int_{a}^{b}f(x)dx =2h\\sum\\limits_{j=0}^{m}f(x_{2j}) + \\frac{h^2(b-a)}{6}f''(\\xi)\\)",
                           style = "font-size:medium"
                         ),
                         
                         br(),
                       ),
                       
                       br(),
                       
                     )
                   ))),
          tabPanel(
            "Ejemplos",
            uiOutput("EICompuestas"),
            mainPanel(navlistPanel(
              tabPanel("Ejemplo 1 ",mainPanel(
                 h4("Fórmulas de integración compuestas", style =
                                   "color: #f2540c"),
                fluidRow(
                  splitLayout(cellWidths = 400,
                             
                            h4(" $$\\int_{0}^{2}|x-2|^3(1-sen(\\pi x))dx $$ "), h4("Considerando los subintervalos: $$[0,\\frac{1}{2}],
                                                                                   [\\frac{1}{2},1],[1,\\frac{3}{2}], [\\frac{3}{2},2]$$")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasIntegracionComp1", "Selecciona un método de integración compuesta:",
                                    c("Regla del trapecio", "Regla de Simpson", "Regla del punto medio"),
                                    selected = "Regla del trapecio",
                                    inline = TRUE
                                  ),
                                  
                                )),
                              
                              
                              
                  ),br()),
                fluidRow(
                  splitLayout(
                    cellWidths = 600,
                    fluidPage( 
                      verticalLayout(
                        uiOutput("DatosIntegracionComp1"),
                        uiOutput("TablaIntegracionComp1"),
                        uiOutput("textoIntegracionComp1"),
                        uiOutput("resultIntegracionComp1"),
                          uiOutput("botonTrapecioComp1"),
                          uiOutput("botonSimpsonComp1"),
                          uiOutput("botonPMComp1")
                        
                      )),
                    plotOutput("graficaIntegracionComp1",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
              tabPanel("Ejemplo 2 ",mainPanel(
                h4("Fórmulas de integración compuestas", style =
                     "color: #f2540c"),
                fluidRow(
                  splitLayout(cellWidths = 400,
                              
                              h4(" $$\\int_{-6}^{0}x^3+6x^2dx $$ "), h4("Considerando los subintervalos: $$[-6,-5],
                                                                                   [-5,-4],[-4,-3],$$$$ [-3,-2],[-2,-1],[-1,0]$$")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasIntegracionComp2", "Selecciona un método de integración compuesta:",
                                    c("Regla del trapecio", "Regla de Simpson", "Regla del punto medio"),
                                    selected = "Regla del trapecio",
                                    inline = TRUE
                                  ),
                                  
                                )),
                              
                              
                              
                  ),br()),
                fluidRow(
                  splitLayout(
                    cellWidths = 600,
                    fluidPage( 
                      verticalLayout(
                        uiOutput("DatosIntegracionComp2"),
                        uiOutput("TablaIntegracionComp2"),
                        uiOutput("textoIntegracionComp2"),
                        uiOutput("resultIntegracionComp2"),
                        uiOutput("AnaliticoIntegracionComp2"),
                        
                        uiOutput("botonTrapecioComp2"),
                        uiOutput("botonSimpsonComp2"),
                        uiOutput("botonPMComp2")
                        
                      )),
                    plotOutput("graficaIntegracionComp2",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
            ))
          )
          
        )
      )
    )
  )
)
v <- reactiveValues(m=data.frame(x = numeric(), "f(x)" = numeric()))
v$m<-NULL
v$m <- data.frame(x = numeric(), "f(x)" = numeric())

inv <- reactiveValues(m=data.frame(x = numeric(), "f(x)" = numeric()))
inv$m<-NULL
inv$m <- data.frame(x = numeric(), "f(x)" = numeric())

her <- reactiveValues(m=data.frame(x = numeric(), "f(x)" = numeric()))
her$m<-NULL
her$m <- data.frame(x = numeric(), "f(x)" = numeric())

server <- function(input, output) {
  ##BIPARTICION
  #EXPONENCIAL
  
  bipartExp <- function(i, n, ab) {
    if (i <= n) {
      output$tablaBipartExp <- renderTable(dfBipartExp[1:i,], digits = 3)
      output$graficaBipartExp <-
        renderPlot(
          gBipartEXP +
            geom_line(
              aes(x = ab[i, 1], colour = "a"),
              color = "#f2540c",
              size = 1
            ) +
            geom_line(
              aes(x = ab[i, 2], colour = "b"),
              color = "#f2540c",
              size = 1
            )
          
          + annotate(
            "text",
            x = ab[i, 1] - 0.05,
            y = -3,
            label = "a",
            size = 5,
            colour = "black",
            vjust = 1.5
          )
          + annotate(
            "text",
            x = ab[i, 2] + 0.05,
            y =  -3,
            label = "b",
            size = 5,
            colour = "black",
            vjust = 1.5
          )
        )
    }
    if (i == n) {
      output$textoBipartExp <-
        renderText(paste(
          withMathJax(
            " La raíz de la ecuación se encuentra entre $$[0.4609375, 0.46484375]$$ con cota de error \\(0.004\\)"
          )
        ))
      output$graficaBipartExp <-
        renderPlot(
          gBipartEXP
          + annotate(
            "point",
            x = 0.4625 ,
            y = 0,
            colour = "red",
            size = 2
          ) +
            annotate(
              "text",
              x = 0.46 + 0.05,
              y = 0,
              label = "Raíz",
              size = 5,
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
  }
  xBipartExp <- seq(0, 1, 0.05)
  dfBipartExp <- bipart_table(function(x)
    (exp(3 * x) - 4), 0, 1, 8, 1)
  colnames (dfBipartExp) <-
    c ('n', 'X_n', 'F (X_n)', '[a, b]', 'Cota de error')
  ab <- bipart(function(x)
    (exp(3 * x) - 4), 0, 1, 8)
  
  xdfBipartExp <- as.data.frame(xBipartExp)
  gBipartEXP <-
    ggplot(xdfBipartExp, aes(x = xBipartExp, y =  exp(3 * xBipartExp) - 4)) +
    annotate(
      "text",
      x = 1,
      y = 12.5,
      label = "f(x)",
      size = 5,
      colour = "black",
      vjust = 1.5
    ) +
    labs(x = "X", y = "F(X)") + theme_bw() +
    geom_line(data = xdfBipartExp) +
    geom_line(color = "blue", size = 1)
  
  
  output$graficaBipartExp <- renderPlot(gBipartEXP)
  
  observeEvent(input$botonBipartExp,
               bipartExp(input$botonBipartExp, 8, ab))
  
  #TRIGONOMETRICA
  
  bipartTrig <- function(i, n, abT) {
    if (i <= n) {
      output$tablaBipartTrig <-
        renderTable(dfBipartTrig[1:i,], digits = 3)
      output$graficaBipartTrig <-
        renderPlot(
          gBipartTRIG +
            geom_line(
              aes(x = abT[i, 1], colour = "a"),
              color = "#f2540c",
              size = 1
            ) +
            geom_line(
              aes(x = abT[i, 2], colour = "b"),
              color = "#f2540c",
              size = 1
            )
          
          + annotate(
            "text",
            x = abT[i, 1] - 0.05 ,
            y = -0.3,
            label = "a",
            size = 5,
            colour = "black",
            vjust = 1.5
          )
          + annotate(
            "text",
            x = abT[i, 2] + 0.05,
            y =  -0.3,
            label = "b",
            size = 5,
            colour = "black",
            vjust = 1.5
          )
        )
    }
    if (i == n) {
      output$textoBipartTrig <-
        renderText(paste(
          withMathJax(
            "  La raíz de la ecuación se encuentra entre $$[-0.8515625, -0.84375]$$ con cota de error \\(0.008\\) "
          )
        ))
      output$graficaBipartTrig <-
        renderPlot(
          gBipartTRIG
          + annotate(
            "point",
            x = -0.849 ,
            y = 0,
            colour = "red",
            size = 2
          ) +
            annotate(
              "text",
              x = -0.849 + 0.1,
              y = 0,
              label = "Raíz",
              size = 5,
              colour = "black",
              vjust = 1.5
            )
        )
    }
  }
  
  xBipartTrig <- seq(-1, 1, 0.05)
  dfBipartTrig <- bipart_table(function(x)
    sin(x) + cos(x ^ 2),-1, 1, 8, 1)
  colnames (dfBipartTrig) <-
    c ('n', 'X_n', 'F (X_n)', '[a, b]', 'Cota de error')
  abT <- bipart(function(x)
    sin(x) + cos(x ^ 2),-1, 1, 8)
  
  
  
  xdfBipartTrig <- as.data.frame(xBipartTrig)
  gBipartTRIG <-
    ggplot(xdfBipartTrig, aes(
      x = xBipartTrig,
      y =  sin(xBipartTrig) + cos(xBipartTrig ^ 2)
    )) +
    labs(x = "X", y = "F(X)") +
    annotate(
      "text",
      x = 0.25 ,
      y = 1.5,
      label = "f(x)",
      colour = "black",
      size = 5,
      vjust = 1.5
    ) + theme_bw() +
    geom_line(data = xdfBipartTrig) +
    geom_line(color = "blue", size = 1)
  
  output$graficaBipartTrig <- renderPlot(gBipartTRIG)
  
  observeEvent(input$botonBipartTrig,
               bipartTrig(input$botonBipartTrig, 8, abT))
  
  
  #PUNTO DE CORTE 2 FUNCIONES
  
  bipartPC <- function(i, n, abPC) {
    if (i <= n) {
      output$tablaBipartPC <- renderTable(dfBipartPC[1:i, ], digits = 3)
      output$graficaBipartPC <-
        renderPlot(
          gBipartPC +
            geom_line(
              aes(x = abPC[i, 1], colour = "a"),
              color = "#f2540c",
              size = 1
            ) +
            geom_line(
              aes(x = abPC[i, 2], colour = "b"),
              color = "#f2540c",
              size = 1
            )
          
          + annotate(
            "text",
            x = abPC[i, 1] - 0.05,
            y = -3,
            label = "a",
            size = 5,
            colour = "black",
            vjust = 1.5
          )
          + annotate(
            "text",
            x = abPC[i, 2] + 0.05,
            y =  -3,
            label = "b",
            size = 5,
            colour = "black",
            vjust = 1.5
          )
        )
    }
    if (i == n) {
      output$textoBipartPC <-
        renderText(paste(
          withMathJax(
            "El punto de corte entre \\(f(x)\\) y \\(g(x)\\) se encuentra en el intervalo $$[3.9921875, 4.01171875]$$ con cota de error \\(0.02\\). En este
                   caso, en $$F(x)=x^2-5x+1=0$$$$x=4 $$ Al evaluar ese valor en $$g(x)=x+1=4+1=5$$ El punto de corte es \\(P(4,5)\\) "
          )
        ))
      output$graficaBipartPC <-
        renderPlot(
          gBipartPC +
            annotate(
              "point",
              x = 4 ,
              y = 5,
              colour = "blue",
              size = 2
            ) +
            annotate(
              "text",
              x = 3.5,
              y = 6.5,
              size = 5,
              label = "Punto de corte",
              colour = "black",
              vjust = 1.5
            )
        )
    }
  }
  
  xBipartPC <- seq(1.5, 7, 0.05)
  dfBipartPC <- bipart_table(function(x)
    (x ^ 2 - 5 * x + 4), 2, 7, 8, 1)
  colnames (dfBipartPC) <-
    c ('n', 'X_n', 'F (X_n)', '[a, b]', 'Cota de error')
  abPC <- bipart(function(x)
    (x ^ 2 - 5 * x + 4), 2, 7, 8)
  
  
  
  xdfBipartPC <- as.data.frame(xBipartPC)
  gBipartPC <-
    ggplot(xdfBipartPC, aes(x = xBipartPC, y =   xBipartPC ^ 2 - 5 * xBipartPC +4)) + 
                        annotate(
                                "text",
                                x = 7,
                                y =  14.5,
                                label = "F(x)",
                                colour = "black",
                                vjust = 1.5,
                                size = 5
                              ) +
    geom_line(aes(x = xBipartPC, y =   (xBipartPC ^ 2 - 4 * xBipartPC + 5)), color = "black", size = 1) + 
                          annotate(
                                "text",
                                x = 7,
                                y =  23,
                                label = "f(x)",
                                colour = "black",
                                vjust = 1.5,
                                size = 5
                          ) +
   geom_line(aes(x = xBipartPC, y =   (xBipartPC + 1)), color = " red", size =  1) + 
                          annotate(
                              "text",
                              x = 7,
                              y =  7,
                              label = "g(x)",
                              colour = "black",
                              vjust = 1.5,
                              size = 5
                            ) +
    labs(x = "X", y = "F(X)") +
    theme_bw() +
    geom_line(data = xdfBipartPC) +
    geom_line(color = "blue", size = 1)
  
  output$graficaBipartPC <- renderPlot(gBipartPC)
  
  observeEvent(input$botonBipartPC, bipartPC(input$botonBipartPC, 8, abPC))
  
  
  #OPTIMIZACION
  
  bipartOpt <- function(i, n, abO) {
    if (i <= n) {
      output$tablaBipartOpt <- renderTable(dfBipartOpt[1:i, ], digits = 3)
      output$graficaBipartOpt <-
        renderPlot(
          gBipartOPT +
            geom_line(
              aes(x = abO[i, 1], colour = "a"),
              color = "#f2540c",
              size = 1
            ) +
            geom_line(
              aes(x = abO[i, 2], colour = "b"),
              color = "#f2540c",
              size = 1
            )
          
          + annotate(
            "text",
            x = abO[i, 1] - 0.05,
            y = -1.5,
            label = "a",
            colour = "black",
            size = 5,
            vjust = 1.5
          )
          + annotate(
            "text",
            x = abO[i, 2] + 0.05,
            y =  -1.5,
            label = "b",
            size = 5,
            colour = "black",
            vjust = 1.5
          )
        )
    }
    if (i == n) {
      output$textoBipartOpt <-
        renderText(
          paste(
            withMathJax(
              " El mínimo o máximo de la función \\(f(x)\\) se encuentra en el intervalo $$[0.5703125, 0.578125]$$ con cota de error \\(0.008\\).En este
                   caso, en $$f'(x)=0$$$$x=0.57736 $$Al evaluar ese valor en $$f''(x)=6x$$$$f''(0.57736)=6·0.57736= 3.46$$"
            ),
            " Se trata de un punto mínimo al ser un valor positivo.<br>",
            withMathJax(
              " El mínimo de la función \\(f(x)\\) es \\(P(0.57736,-1.3849)\\) "
            )
          )
        )
      output$graficaBipartOpt <-
        renderPlot(
          gBipartOPT +
            annotate(
              "point",
              x = 0.57736 ,
              y = -1.3849,
              colour = "red",
              size = 2
            ) +
            annotate(
              "text",
              x = 0.57736 + 0.1,
              y = -1.5,
              label = "Punto mínimo",
              size = 5,
              colour = "black",
              vjust = 1.5
            )
        )
    }
  }
  
  xBipartOpt <- seq(0, 2, 0.05)
  dfBipartOpt <- bipart_table(function(x)
    (3 * x ^ 2 - 1), 0, 2, 8, 1)
  colnames (dfBipartOpt) <-
    c ('n', 'X_n', 'F (X_n)', '[a, b]', 'Cota de error')
  abO <- bipart(function(x)
    (3 * x ^ 2 - 1), 0, 2, 8)
  
  xdfBipartOpt <- as.data.frame(xBipartOpt)
  gBipartOPT <-
    ggplot(xdfBipartOpt, aes(x = xBipartOpt, y =   (3 * xBipartOpt ^ 2 - 1))) + annotate(
      "text",
      x = 2,
      y =  4.5,
      label = "f(x)",
      colour = "black",
      vjust = 1.5,
      size = 5
    ) + geom_line(aes(x = xBipartOpt, y =   (xBipartOpt ^ 3 - xBipartOpt - 1)), color =
                    "black", size = 1) + annotate(
                      "text",
                      x = 2,
                      y =  10,
                      label = "f '(x)",
                      colour = "black",
                      vjust = 1.5,
                      size = 5
                    ) +
    labs(x = "X", y = "F(X)'") +
    theme_bw() +
    geom_line(data = xdfBipartOpt) +
    geom_line(color = "blue", size = 1)
  
  output$graficaBipartOpt <- renderPlot(gBipartOPT)
  observeEvent(input$botonBipartOpt,
               bipartOpt(input$botonBipartOpt, 8, abO))
  
  
  
  ##NEWTON
  #EXPONENCIAL
  newtonExp <- function(i, n) {
    if (i <= n) {
      output$tablaNewtonExp <- renderTable(dfNewtonExp2[1:i, ], digits = 3)
      output$textoNewtonExp <-
        renderText(paste(
          withMathJax(
            "$$x_",
            i,
            "= ",
            round(dfNewtonExp2[i, 2], 2),
            "-\\frac{",
            round(dfNewtonExp2[i, 3], 2) ,
            "}{",
            round(dfNewtonExp2[i, 4], 2),
            "} = ",
            round(dfNewtonExp2[i, 5], 2) ,
            "$$"
          )
        ))
      output$graficaNewtonExp <-
        renderPlot(
          g +   geom_point(aes(x = round(
            dfNewtonExp[i + 1, 2], 2
          )), color = "black", size = 1) +
            
            geom_abline(
              slope = slope(round(dfNewtonExp[i, 2], 2)),
              intercept = intercept(round(dfNewtonExp[i, 2], 2)),
              color = "#f2540c",
              size = 1
            )
          + annotate(
            "point",
            x = round(dfNewtonExp[i, 2], 2) ,
            y = round(dfNewtonExp[i, 4], 2),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonExp[i, 2], 2) - 0.1,
              y = round(dfNewtonExp[i, 4], 2),
              label = round(dfNewtonExp[i, 2], 2),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfNewtonExp[i + 1, 2], 2) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonExp[i + 1, 2], 2) + 0.15,
              y = 0,
              label = round(dfNewtonExp[i + 1, 2], 2),
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    if (i == n) {
      output$textoNewtonExp <-
        renderText(paste(withMathJax(
          "La raíz de la ecuación es \\(1.058\\) "
        )))
      output$graficaNewtonExp <-
        renderPlot(
          g
          + annotate(
            "point",
            x = 1.058006 ,
            y = 0,
            colour = "red",
            size = 2
          ) +
            annotate(
              "text",
              x = 1.058006 + 0.1,
              y = 0,
              label = "Raíz",
              size = 5,
              colour = "black",
              vjust = 1.5
            )
        )
    }
    
    
  }
  
  dfNewtonExp <- newton(function(x)
    (exp(x) + x ^ 2 - 4),  function(x)
      (exp(x) + 2 * x)
    , 3.5)
  dfNewtonExp2 <- newtonII(function(x)
    (exp(x) + x ^ 2 - 4),  function(x)
      (exp(x) + 2 * x)
    , 3.5)
  colnames (dfNewtonExp2) <-
    c ('k', 'X_k', 'F(X_k)', 'F`(X_k)', 'X_(k+1)')
  
  
  #GRAFICA
  xNewtonExp <- seq(0, 3.5, 0.05)
  xdfNewtonExp <- as.data.frame(xNewtonExp)
  g <-
    ggplot(xdfNewtonExp, aes(x = xNewtonExp, y =  exp(xNewtonExp) + xNewtonExp ^ 2 - 4)) +
    labs(x = "X", y = "F(X)") +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "f(x)",
      colour = "black",
      vjust = 1.5,
      size = 5
    ) + theme_bw() +
    geom_line(data = xdfNewtonExp) +
    geom_line(color = "blue", size = 1)
  
  
  #Y=mx+b  m=> SLOPE b=>INTERCEPT
  slope <- function(x) {
    xincrement <- 0.2
    yincrement <-
      (exp((x + 0.1)) + (x + 0.1) ^ 2 - 4) - (exp((x - 0.1)) + (x - 0.1) ^ 2 - 4)
    slope <- yincrement / xincrement
    return(slope)
  }
  
  intercept <- function(x) {
    intercept <- (exp(x) + x ^ 2 - 4) - slope(x) * x
    yvalues <- slope(x) * xNewtonExp + intercept
    xdfNewtonExp <- as.data.frame(cbind(xdfNewtonExp, yvalues))
    return(intercept)
  }
  
  output$graficaNewtonExp <- renderPlot(g)
  
  observeEvent(input$botonNewtonExp, newtonExp(input$botonNewtonExp, 6))
  
  
  #TRIGONOMÉTRICA
  newtonTrig <- function(i, n) {
    if (i <= n) {
      output$tablaNewtonTrig <- renderTable(dfNewtonTrig2[1:i, ], digits = 4)
      output$textoNewtonTrig <-
        renderText(paste(
          withMathJax(
            "$$x_",
            i,
            "= ",
            round(dfNewtonTrig2[i, 2], 4),
            "-\\frac{",
            round(dfNewtonTrig2[i, 3], 2) ,
            "}{",
            round(dfNewtonTrig2[i, 4], 2),
            "} = ",
            round(dfNewtonTrig2[i, 5], 2) ,
            "$$"
          )
        ))
      output$graficaNewtonTrig <-
        renderPlot(
          gTrig + geom_point(aes(x = round(
            dfNewtonTrig[i + 1, 2], 2
          )), color = "black", size = 1) +
            geom_abline(
              slope = slopeTrig(dfNewtonTrig[i, 2]),
              intercept = interceptTrig(dfNewtonTrig[i, 2]),
              color = "#f2540c",
              size = 1
            )
          + annotate(
            "point",
            x = round(dfNewtonTrig[i, 2], 4) ,
            y = round(dfNewtonTrig[i, 4], 4),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonTrig[i, 2], 4) - 0.1,
              y = round(dfNewtonTrig[i, 4], 4),
              label = round(dfNewtonTrig[i, 2], 4),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfNewtonTrig[i + 1, 2], 4) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonTrig[i + 1, 2], 4) + 0.15,
              y = 0,
              label = round(dfNewtonTrig[i + 1, 2], 4),
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    if (i == n) {
      output$textoNewtonTrig <-
        renderText(paste(withMathJax(
          " La raíz de la ecuación es \\(0.824\\)"
        )))
      output$graficaNewtonTrig <-
        renderPlot(
          gTrig
          + annotate(
            "point",
            x = 0.824132 ,
            y = 0,
            colour = "red",
            size = 2
          ) +
            annotate(
              "text",
              x = 0.824132 + 0.1,
              y = 0,
              label = "Raíz",
              size = 5,
              colour = "black",
              vjust = 1.5
            )
        )
    }
    
  }
  dfNewtonTrig <- newton(function(x)
    (x ^ 2 - cos(x)),  function(x)
      (2 * x + sin(x))
    , 2)
  dfNewtonTrig2 <- newtonII(function(x)
    (x ^ 2 - cos(x)),  function(x)
      (2 * x + sin(x))
    , 2)
  colnames (dfNewtonTrig2) <-
    c ('k', 'X_k', 'F(X_k)', 'F`(X_k)', 'X_(k+1)')
  
  #GRAFICA
  xNewtonTrig <- seq(0, 2, 0.05)
  xdfNewtonTrig <- as.data.frame(xNewtonTrig)
  gTrig <-
    ggplot(xdfNewtonTrig, aes(x = xNewtonTrig, y = (xNewtonTrig ^ 2 - cos(xNewtonTrig)))) +
    labs(x = "X", y = "F(X)") +
    annotate(
      "text",
      x = 0,
      y = -0.5,
      label = "f(x)",
      size = 5,
      colour = "black",
      vjust = 1.5
    ) +
    theme_bw() +
    geom_line(data = xdfNewtonTrig) +
    geom_line(color = "blue", size = 1)
  
  
  #Y=mx+b  m=> SLOPE b=>INTERCEPT
  slopeTrig <- function(x) {
    xincrementT <- 0.2
    yincrementT <- ((x + 0.1) ^ 2 - cos(x + 0.1)) - ((x - 0.1) ^ 2 - cos(x -
                                                                           0.1))
    slopeT <- yincrementT / xincrementT
    return(slopeT)
  }
  
  interceptTrig <- function(x) {
    interceptT <- (x ^ 2 - cos(x)) - slopeTrig(x) * x
    yvaluesT <- slopeTrig(x) * xNewtonTrig + interceptT
    xdfNewtonTrig <- as.data.frame(cbind(xdfNewtonTrig, yvaluesT))
    return(interceptT)
  }
  
  output$graficaNewtonTrig <- renderPlot(gTrig)
  
  observeEvent(input$botonNewtonTrig,
               newtonTrig(input$botonNewtonTrig, 5))
  
  #PUNTO DE CORTE
  newtonPC <- function(i, n) {
    if (i < 4) {
      output$tablaNewtonPC <- renderTable(dfNewtonPC2[1:i, ], digits = 4)
      output$textoNewtonPC <-
        renderText(paste(
          withMathJax(
            "$$x_",
            i,
            "= ",
            round(dfNewtonPC2[i, 2], 4),
            "-\\frac{",
            round(dfNewtonPC2[i, 3], 2) ,
            "}{",
            round(dfNewtonPC2[i, 4], 2),
            "} = ",
            round(dfNewtonPC2[i, 5], 2) ,
            "$$"
          )
        ))
      output$graficaNewtonPC <-
        renderPlot(
          gPC + geom_point(aes(x = round(
            dfNewtonPC[i + 1, 2], 2
          )), color = "black", size = 1) +
            geom_abline(
              slope = slopePC(dfNewtonPC[i, 2]),
              intercept = interceptPC(dfNewtonPC[i, 2]),
              color = "#f2540c",
              size = 1
            )
          + annotate(
            "point",
            x = round(dfNewtonPC[i, 2], 4) ,
            y = round(dfNewtonPC[i, 4], 4),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonPC[i, 2], 4) + 0.2,
              y = round(dfNewtonPC[i, 4], 4) + 0.3,
              label = round(dfNewtonPC[i, 2], 4),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfNewtonPC[i + 1, 2], 4) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonPC[i + 1, 2], 4) - 0.2,
              y = -0.12,
              label = round(dfNewtonPC[i + 1, 2], 4),
              colour = "black",
              vjust = 1.5
            )
        )
      
    } else if (i > 3 && i < n) {
      output$tablaNewtonPC <- renderTable(dfNewtonPC2[1:i, ], digits = 4)
      output$textoNewtonPC <-
        renderText(paste(
          withMathJax(
            "$$ x_",
            i,
            "= ",
            round(dfNewtonPC2[i, 2], 4),
            "-\\frac{",
            round(dfNewtonPC2[i, 3], 2) ,
            "}{",
            round(dfNewtonPC2[i, 4], 2),
            "} = ",
            round(dfNewtonPC2[i, 5], 2) ,
            "$$"
          )
        ))
      output$graficaNewtonPC <-
        renderPlot(
          gPC + geom_point(aes(x = round(
            dfNewtonPC[i + 1, 2], 2
          )), color = "black", size = 1) +
            geom_abline(
              slope = slopePC(dfNewtonPC[i, 2]),
              intercept = interceptPC(dfNewtonPC[i, 2]),
              color = "#f2540c",
              size = 1
            )
          + annotate(
            "point",
            x = round(dfNewtonPC[i, 2], 4) ,
            y = round(dfNewtonPC[i, 4], 4),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonPC[i, 2], 4) - 0.5,
              y = round(dfNewtonPC[i, 4], 4) - 0.7,
              label = round(dfNewtonPC[i, 2], 4),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfNewtonPC[i + 1, 2], 4) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonPC[i + 1, 2], 4) + 0.4,
              y = 3,
              label = round(dfNewtonPC[i + 1, 2], 4),
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    if (i == n) {
      output$textoNewtonPC <-
        renderText(paste(
          withMathJax("Hemos calculado que \\(x = -1.5348\\)"),
          "<br>",
          withMathJax(
            " Evaluamos ese valor en la función:$$g(x)=sen(-1.5348)-(-1.5348)^3=2.616$$El punto de corte entre \\(f(x)\\) y \\(g(x)\\) es \\(P(-1.5348,2.616)\\)"
          )
        ))
      output$graficaNewtonPC <-
        renderPlot(
          gPC  + geom_line(
            aes(
              x = xNewtonPC,
              y =   (xNewtonPC ^ 4 + 3 * xNewtonPC ^ 2 - 10)
            ),
            color = "black",
            size = 1
          ) + annotate(
            "text",
            x = -3.5,
            y =  150,
            label = "f(x)",
            colour = "black",
            vjust = 1.5,
            size = 5
          ) + geom_line(
            aes(x = xNewtonPC, y =   (sin(xNewtonPC) - xNewtonPC ^ 3)),
            color = " red",
            size = 1
          ) + annotate(
            "text",
            x = -3.5,
            y =  30,
            label = "g(x)",
            colour = "black",
            vjust = 1.5,
            size = 5
          )
          
          + annotate(
            "point",
            x = -1.5348 ,
            y = 2.616,
            colour = "blue",
            size = 2
          ) +
            annotate(
              "text",
              x = -1.5348 + 0.5,
              y = 10,
              label = "Punto de corte",
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    
    
  }
  
  
  dfNewtonPC <- newton(function(x)
    (x ^ 4 + x ^ 3 + 3 * x ^ 2 - sin(x) - 10),  function(x)
      (4 * x ^ 3 + 3 * x ^ 2 + 6 * x - cos(x))
    ,-0.5)
  dfNewtonPC2 <- newtonII(function(x)
    (x ^ 4 + x ^ 3 + 3 * x ^ 2 - sin(x) - 10),  function(x)
      (4 * x ^ 3 + 3 * x ^ 2 + 6 * x - cos(x))
    ,-0.5)
  colnames (dfNewtonPC2) <-
    c ('k', 'X_k', 'F(X_k)', 'F`(X_k)', 'X_(k+1)')
  
  #GRAFICA
  xNewtonPC <- seq(-3.5, 1, 0.05)
  xdfNewtonPC <- as.data.frame(xNewtonPC)
  gPC <-
    ggplot(xdfNewtonPC,
           aes(
             x = xNewtonPC,
             y = xNewtonPC ^ 4 + xNewtonPC ^ 3 + 3 * xNewtonPC ^ 2 - sin(xNewtonPC) -
               10
           )) + annotate(
             "text",
             x = -3.5,
             y =  110,
             label = "F(x)",
             colour = "black",
             vjust = 1.5,
             size = 5
           ) +
    labs(x = "X", y = "F(X)") +
    theme_bw() +
    geom_line(data = xdfNewtonPC) +
    geom_line(color = "blue", size = 1)
  
  
  #Y=mx+b  m=> SLOPE b=>INTERCEPT
  slopePC <- function(x) {
    xincrementPC <- 0.2
    yincrementPC <-
      ((x + 0.1) ^ 4 + (x + 0.1) ^ 3 + 3 * (x + 0.1) ^ 2 - sin((x + 0.1)) -
         10) - ((x - 0.1) ^ 4 + (x - 0.1) ^ 3 + 3 * (x - 0.1) ^ 2 - sin((x - 0.1)) -
                  10)
    slopepto <- yincrementPC / xincrementPC
    return(slopepto)
  }
  
  interceptPC <- function(x) {
    interceptpto <- (x ^ 4 + x ^ 3 + 3 * x ^ 2 - sin(x) - 10) - slopePC(x) * x
    yvaluesPC <- slopePC(x) * xNewtonPC + interceptpto
    xdfNewtonPC <- as.data.frame(cbind(xdfNewtonPC, yvaluesPC))
    return(interceptpto)
  }
  
  output$graficaNewtonPC <- renderPlot(gPC)
  
  observeEvent(input$botonNewtonPC, newtonPC(input$botonNewtonPC, 7))
  
  
  
  #OPTIMIZACIÓN
  newtonOpt <- function(i, n) {
    if (i <= n) {
      output$tablaNewtonOpt <- renderTable(dfNewtonOpt2[1:i, ], digits = 4)
      output$textoNewtonOpt <-
        renderText(paste(
          withMathJax(
            "$$x_{",
            i,
            "}= ",
            round(dfNewtonOpt2[i, 2], 4),
            "-\\frac{",
            round(dfNewtonOpt2[i, 3], 2) ,
            "}{",
            round(dfNewtonOpt2[i, 4], 2),
            "} = ",
            round(dfNewtonOpt2[i, 5], 2),
            "$$"
          )
        ))
      output$graficaNewtonOpt <-
        renderPlot(
          gOpt + geom_point(aes(x = round(
            dfNewtonOpt[i + 1, 2], 2
          )), color = "black", size = 1) +
            geom_abline(
              slope = slopeOpt(dfNewtonOpt[i, 2]),
              intercept = interceptOpt(dfNewtonOpt[i, 2]),
              color = "#f2540c",
              size = 1
            )
          + annotate(
            "point",
            x = round(dfNewtonOpt[i, 2], 4) ,
            y = round(dfNewtonOpt[i, 4], 4),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonOpt[i, 2], 4) - 0.1,
              y = round(dfNewtonOpt[i, 4], 4),
              label = round(dfNewtonOpt[i, 2], 4),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfNewtonOpt[i + 1, 2], 4) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonOpt[i + 1, 2], 4) + 0.15,
              y = 0,
              label = round(dfNewtonOpt[i + 1, 2], 4),
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    if (i == n) {
      output$textoNewtonOpt <-
        renderText(
          paste(
            withMathJax("$$x_5=0 $$$$  f''(0)=6·0+8=8 $$"),
            "El resultado positivo implica un punto mínimo <br>",
            withMathJax("El mínimo de la funcion \\(f(x)\\) es \\(P(0,-10)\\)")
          )
        )
      output$graficaNewtonOpt <-
        renderPlot(
          gOpt  + geom_abline(
            slope = slopeOpt(dfNewtonOpt[i, 2]),
            intercept = interceptOpt(dfNewtonOpt[i, 2]),
            color = "#f2540c",
            size = 1
          )
          + annotate(
            "point",
            x = round(dfNewtonOpt[i, 2], 4) ,
            y = round(dfNewtonOpt[i, 4], 4),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonOpt[i, 2], 4) - 0.1,
              y = round(dfNewtonOpt[i, 4], 4),
              label = round(dfNewtonOpt[i, 2], 4),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfNewtonOpt[i + 1, 2], 4) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfNewtonOpt[i + 1, 2], 4) + 0.15,
              y = 0,
              label = round(dfNewtonOpt[i + 1, 2], 4),
              colour = "black",
              vjust = 1.5
            )
          + annotate(
            "point",
            x = 0 ,
            y = -10,
            colour = "red",
            size = 2
          ) +
            annotate(
              "text",
              x = 0 + 0.15,
              y = -10,
              label = "Punto mínimo",
              size = 5,
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    
    
  }
  
  dfNewtonOpt <- newton(function(x)
    (3 * x ^ 2 + 8 * x),  function(x)
      (6 * x + 8)
    , 0.75)
  dfNewtonOpt2 <- newtonII(function(x)
    (3 * x ^ 2 + 8 * x),  function(x)
      (6 * x + 8)
    , 0.75)
  colnames (dfNewtonOpt2) <-
    c ('k', 'X_k', 'F(X_k)', 'F``(X_k)', 'X_(k+1)')
  
  xNewtonOpt <- seq(-1.5, 1.5, 0.05)
  
  
  #GRAFICA
  
  xNewtonOptdf <- as.data.frame(xNewtonOpt)
  gOpt <-
    ggplot(xNewtonOptdf, aes(x = xNewtonOpt, y = (3 * xNewtonOpt ^ 2 + 8 *(xNewtonOpt)))) +
                            annotate(
                              "text",
                              x = 1.5,
                              y =  16,
                              label = "f '(x)",
                              colour = "black",
                              vjust = 1.5,
                              size = 5
                          ) +
    geom_line(aes(x = xNewtonOpt, y =  (xNewtonOpt ^ 3 + 4 * (xNewtonOpt) ^2 - 10)), color = "black", size = 1) + 
                          annotate(
                            "text",
                            x = 1.5,
                            y =  -2,
                            label = "f(x)",
                            colour = "black",
                            vjust = 1.5,
                            size = 5
                          ) +
    labs(x = "X", y = "F(X)") +
    theme_bw() +
    geom_line(data = xNewtonOptdf) +
    geom_line(color = "blue", size = 1)
  
  
  #Y=mx+b  m=> SLOPE b=>INTERCEPT
  slopeOpt <- function(x) {
    xincrementO <- 0.2
    yincrementO <-
      (3 * (x + 0.1) ^ 2 + 8 * (x + 0.1)) - (3 * (x - 0.1) ^ 2 + 8 * (x -
                                                                        0.1))
    slopeO <- yincrementO / xincrementO
    return(slopeO)
  }
  
  interceptOpt <- function(x) {
    interceptO <- (3 * x ^ 2 + 8 * (x)) - slopeOpt(x) * x
    yvaluesO <- slopeOpt(x) * xNewtonOpt + interceptO
    xNewtonOptdf <- as.data.frame(cbind(xNewtonOptdf, yvaluesO))
    return(interceptO)
  }
  
  output$graficaNewtonOpt <- renderPlot(gOpt)
  
  observeEvent(input$botonNewtonOpt, newtonOpt(input$botonNewtonOpt, 5))
  
  
  
  ##INTERPOLACION NEWTON----------------------------------------------------------------------------------
  
  ##INTRODUCIR DATOS
  observe( v$m <- data.frame(x = numeric(), "f(x)" = numeric()))
  cont<-0
  observeEvent(input$botonNewtonPuntos,{
    req(input$x,input$y) 
    
    tmp <- data.frame(x = input$x,y  = input$y)
    
    colnames (tmp) <- c ( 'x', 'f(x)')
    
    size<-dim(v$m)
    if(size[2]<3){
      output$textoInterpolacionNewton<-renderText(paste(""))
      
      v$m <- rbind(v$m,tmp)
      output$InterpolacionNewton <- DT::renderDataTable({
        v$m
      },options = list(dom = 't')) 
      gInterN<-ggplot((v$m)) +
        geom_point(aes(x =v$m[,1],y=v$m[,2]))+theme_bw()+ labs(x = "X", y = "F(X)")
      output$graficaInterpolacionNewton <-
        renderPlot(gInterN)
    }
  })
  
  ##REALIZAR TABLA
  
  
  counter <- reactiveValues(countervalue = 0)
  
  ##RESET
  observeEvent(input$botonNewtonReset,{
  v$m <- data.frame(x = numeric(), "f(x)" = numeric())
  gInterInf<-ggplot( v$m )
  output$graficaInterpolacionNewton <-NULL
  output$textoInterpolacionNewton<-NULL
  counter$countervalue <-0 
  })
  
  ##ELIMINAR DATOS QUE SOBRAN
  observeEvent(input$botonInterpolacionNewton,{
    tryCatch(
      expr = {
    counter$countervalue <- counter$countervalue + 1   
    tabla<- as.data.frame(v$m)
    
    tam<-dim(tabla)
    
    if(tam[2]<tam[1]+1){
      
      tab<- intNewton_table(tabla[,1],tabla[,2])
      tam2<-dim(tab)  
      bool<-(tab[,1]== pracma::zeros(tam2[1],1))
      cont<-0
      for(j in 1:tam2[1]){
        if(bool[j,1]==TRUE){
          cont<-cont+1
        
        }
      }
      
     if(cont==tam2[1]){ 
        v$m <- data.frame(x = numeric(), "f(x)" = numeric())
        gInterInf<-ggplot( v$m )
        output$graficaInterpolacionNewton <-NULL
        output$textoInterpolacionNewton<-renderText(paste("<h4 style=color:#FF0000;>ERROR <h4> <h4>Datos no validos, por favor introduce nuevos puntos. <h4>"))
        counter$countervalue <-0  
      }else{
        fx_n<- tab[, counter$countervalue+2]
        v$m <- cbind(v$m,
                     fx_n)
        colnames(v$m)[counter$countervalue+2]<-paste(counter$countervalue,"º dif. divididas")
        
        output$InterpolacionNewton <- DT::renderDataTable({
          v$m
        },options = list(dom = 't'))
        
      }
    }else{
      polinomio<-polinomioInt(v$m[1,],v$m[,1])
      output$textoInterpolacionNewton<-renderText(paste(withMathJax("El polinomio de interpolación es: $$P(x) =", polinomio,"$$" )))
    }
      },
    error=function(e){
      v$m <- data.frame(x = numeric(), "f(x)" = numeric())
      gInterInf<-ggplot( v$m )
      output$graficaInterpolacionNewton <-NULL
      output$textoInterpolacionNewton<-renderText(paste("<h4 style=color:#FF0000;>ERROR <h4> <h4>Datos no validos, por favor introduce nuevos puntos. <h4>"))
      counter$countervalue <-0  
      
    }
    )
    
  })
  
  
  ##INTERPOLACION INVERSA
  
  observe( inv$m <- data.frame(x = numeric(), "f(x)" = numeric()))

  observeEvent(input$botonInversaPuntos,{
    
    req(input$x1,input$y1) 
    
    tmpInv <- data.frame(x = input$x1,y  = input$y1)
    colnames (tmpInv) <- c ( 'x', 'f(x)')
    
    sizeInv<-dim(inv$m)
    if(sizeInv[2]<3){
      output$textoInterpolacionInversa<-renderText(paste(""))
      
      inv$m <- rbind(inv$m,tmpInv)
      output$InterpolacionInversa <- DT::renderDataTable({
        inv$m
      },options = list(dom = 't')) 
      gInterInv<-ggplot((inv$m)) +
        geom_point(aes(x =inv$m[,1],y=inv$m[,2]))+theme_bw()+ labs(x = "X", y = "F(X)")
      output$graficaInterpolacionInversa <-
        renderPlot(gInterInv)
    }
  })
  
  ##REALIZAR TABLA
  
  
  counterInv <- reactiveValues(countervalue = 0,countervalueInv2=0)
  
  ##RESET
  observeEvent(input$botonInversaReset,{
    inv$m <- data.frame(x = numeric(), "f(x)" = numeric())
    gInterInf2<-ggplot( inv$m )
    output$graficaInterpolacionInversa <-NULL
    output$textoInterpolacionInversa<-NULL
    output$textoInterpolacionValor<-NULL
    output$y2<-NULL
    output$botonIntInvV<-NULL
    counterInv$countervalue <-0 
    counterInv$countervalueinv2 <-0 
  })
  
  ##ELIMINAR DATOS QUE SOBRAN
  observeEvent(input$botonInterpolacionInversa, {
    counterInv$countervalue <- counterInv$countervalue + 1
    tryCatch(
      expr = {
        tablaInv <- as.data.frame(inv$m)
        tamInv <- dim(tablaInv)
        
        if (tamInv[2] < tamInv[1] + 1) {
          tabInv <- intNewton_table(tablaInv[, 1], tablaInv[, 2])
          tamInv2 <- dim(tabInv)
          
          bool2 <- (tabInv[, 1] == pracma::zeros(tamInv2[1], 1))
          cont2 <- 0
          for (j in 1:tamInv2[1]) {
            if (bool2[j, 1] == TRUE) {
              cont2 <- cont2 + 1
              
            }
          }
          
          if (cont2 == tamInv2[1]) {
            inv$m <- data.frame(x = numeric(), "f(x)" = numeric())
            gInterInf2 <- ggplot(inv$m)
            output$graficaInterpolacionInversa <- NULL
            output$textoInterpolacionInversa <-
              renderText(
                paste(
                  "<h4 style=color:#FF0000;>ERROR <h4> <h4>Datos no validos, por favor introduce nuevos puntos. <h4>"
                )
              )
            counterInv$countervalue <- 0
          } else{
            fx_n2 <- tabInv[, counterInv$countervalue + 2]
            inv$m <- cbind(inv$m,
                           fx_n2)
            colnames(inv$m)[counterInv$countervalue + 2] <-
              paste(counterInv$countervalue ,"º dif. divididas")
            
            output$InterpolacionInversa <- DT::renderDataTable({
              inv$m
            }, options = list(dom = 't'))
            
          }
        } else{
          polinomioInv <- polinomioInt(inv$m[1, ], inv$m[, 1])
          
          output$textoInterpolacionInversa <-
            renderText(paste(
              withMathJax(
                "El polinomio de interpolación es \\(P(x) =",
                polinomioInv,
                "\\)"
              )
            ))
          output$textoInterpolacionValor <-
            renderText(paste(
              withMathJax("Introduce el valor para resolver la ecuación : ")
            ))
          output$y2 <- renderUI({
            numericInput("y2", "", 0, min = -50, max = 50)
          })
          output$botonIntInvV <- renderUI({
            actionButton("botonIntInvV", "Introducir valor", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
          })
        }
      },
      error = function(e) {
        inv$m <- data.frame(x = numeric(), "f(x)" = numeric())
        gInterInf2 <- ggplot(inv$m)
        output$graficaInterpolacionInversa <- NULL
        output$textoInterpolacionInversa <-
          renderText(
            paste(
              "<h4 style=color:#FF0000;>ERROR <h4> <h4>Valor no válido o polinomio sin solución, por favor introduce nuevo valor o nuevos puntos. <h4>"
            )
          )
        counterInv$countervalue <- 0
      }
    )
    
    
  })
  
  observeEvent(input$botonIntInvV,{
    counterInv$countervalueInv2 <- counterInv$countervalueInv2 + 1   
    req(input$y2)
    y2<-input$y2
   
      tablaInv2<- as.data.frame(inv$m )
      tamInv2<-dim(tablaInv2)
      polinomioInv2<-polinomioInt(inv$m [1,],inv$m [,1])
      if(y2<0){
      poly <- ysym(paste(polinomioInt2(inv$m[1,],inv$m[,1]),"+",y2))
      }else{
        poly <- ysym(paste(polinomioInt2(inv$m[1,],inv$m[,1]),"-",y2))
    } 
      print(poly)
      tryCatch(
        expr = {
          resultado<-solve(poly, "x")
          print(resultado)
          output$textoInterpolacionInversa<-renderText(paste(withMathJax("El polinomio de interpolación es \\(P(x) =", polinomioInv2,"\\) $$$$ Resolvemos la ecuación : \\(",polinomioInv2,"=",y2,"\\) $$$$ El resultado es \\(", as.character(resultado),"\\)")))
        },
        error = function(e){ 
          inv$m <- data.frame(x = numeric(), "f(x)" = numeric())
          gInterInf2<-ggplot( inv$m )
          output$graficaInterpolacionInversa <-NULL
          output$textoInterpolacionInversa<-renderText(paste("<h4 style=color:#FF0000;>ERROR <h4> <h4>Datos no validos, por favor introduce nuevos puntos. <h4>"))
          output$y2<-NULL
          output$botonIntInvV<-NULL
          counterInv$countervalue <-0 
          counterInv$countervalueinv2 <-0
          output$textoInterpolacionValor<-NULL
          
        }
       
      )
       
    
  })
 
  #OSCULATORIA HERMITE
  
  
  ##INTRODUCIR DATOS
  observe( her$m <- data.frame(x = numeric(), "f(x)" = numeric()))
  observeEvent(input$botonHermitePuntos,{
  
    req(input$xh,input$yh,input$yhd) 
   
    tmpH<- data.frame(x = input$xh,y  = input$yh, y2=input$yhd)
    
    colnames (tmpH) <- c ( 'x', 'f(x)','f(x)`')
    
    sizeH<-dim(her$m)
    if(sizeH[2]<4){
      output$textoInterpolacionHermite<-renderText(paste(""))
      
      her$m <- rbind(her$m,tmpH)
      
      output$InterpolacionHermite <- DT::renderDataTable({
        her$m
      },options = list(dom = 't')) 
      gInterN<-ggplot((her$m)) +
        geom_point(aes(x =her$m[,1],y=her$m[,2]))+theme_bw()+ labs(x = "X", y = "F(X)")
      output$graficaInterpolacionHermite <-
        renderPlot(gInterN)
    }
  })
  
  ##REALIZAR TABLA
  
  
  counterH <- reactiveValues(countervalue = 0)
  
  ##RESET
  observeEvent(input$botonHermiteReset,{
    her$m <- data.frame(x = numeric(), "f(x)" = numeric())
    gInterInf<-ggplot( her$m )
    output$graficaInterpolacionHermite <-NULL
    output$textoInterpolacionHermite<-NULL
    counterH$countervalue <-0 
  })
  
  ##ELIMINAR DATOS QUE SOBRAN
  observeEvent(input$botonInterpolacionHermite,{
    tryCatch(
      expr = {
        counterH$countervalue <- counterH$countervalue + 1   
        tablaH<- as.data.frame(her$m)
        print(tablaH)
        tamH<-dim(tablaH)
        print(tamH)
        if(tamH[2]<tamH[1]+2){
          print("hola")
          tabH<- hermite_table(tablaH[,1],tablaH[,2],tablaH[,3],3)
          print(tabH)
          tamH2<-dim(tabH)  
          boolH<-(tabH[,1]== pracma::zeros(tamH2[1],1))
          contH<-0
          for(j in 1:tamH2[1]){
            if(boolH[j,1]==TRUE){
              contH<-contH+1
              
            }
          }
          
          if(contH==tamH2[1]){ 
            her$m <- data.frame(x = numeric(), "f(x)" = numeric())
            gInterInfH<-ggplot( her$m )
            output$graficaInterpolacionHermite <-NULL
            output$textoInterpolacionHermite<-renderText(paste("<h4 style=color:#FF0000;>ERROR <h4> <h4>Datos no validos, por favor introduce nuevos puntos. <h4>"))
            counterH$countervalue <-0  
          }else{
            fx_h<- tab[, counterH$countervalue+2]
            her$m <- cbind(her$m,
                         fx_h)
            colnames(her$m)[counterH$countervalue+2]<-paste(counterH$countervalue,"º dif. divididas")
            
            output$InterpolacionHermite <- DT::renderDataTable({
              her$m
            },options = list(dom = 't'))
            
          }
        }else{
          polinomioH<-polinomioInt(her$m[1,],her$m[,1])
          output$textoInterpolacionHermite<-renderText(paste(withMathJax("El polinomio de interpolación es: $$P(x) =", polinomioH,"$$" )))
        }
      },
      error=function(e){
        her$m <- data.frame(x = numeric(), "f(x)" = numeric())
        gInterInfH<-ggplot( her$m )
        output$graficaInterpolacionHermite <-NULL
        output$textoInterpolacionHermite<-renderText(paste("<h4 style=color:#FF0000;>ERROR <h4> <h4>Datos no validos, por favor introduce nuevos puntos. <h4>"))
        counterH$countervalue <-0  
        
      }
    )
    
  })
  
  
  
  
  ##DERIVACIÓN COMPLICADA
  observeEvent(ignoreInit = FALSE, c(
   input$metodosDerivacion1, input$h
  ),{
    req(input$h)
    
    derC<- derivC(function(x)
      x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3,input$h)
    derD<- derivD(function(x)
      x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3,input$h)
    derI<- derivI(function(x)
      x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3,input$h)
    h<-input$h
    
    
    if(input$metodosDerivacion1== "Centrada" ){
      output$textoDerivacion1 <- renderText( paste( withMathJax("$$f'(1.3) \\approx \\frac{f(1.3+",h,")-f(1.3-",h,")}{2·",h,"} \\approx$$$$\\frac{ f(",1.3+h,")-f(",1.3-h,")}{",2*h,"}\\approx \\frac{ ",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3+h),2),
                                                        "-",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3-h),2),"}{",2*h,"}\\approx",derC,"$$")))
      
    }else if(input$metodosDerivacion1== "Descentrada izquierda" ){
      output$textoDerivacion1 <- renderText( paste( withMathJax("$$f'(1.3) \\approx \\frac{f(1.3)-f(1.3-",h,")}{",h,"} \\approx$$$$\\frac{ f(",1.3,")-f(",1.3-h,")}{",h,"}\\approx \\frac{ ",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3),2),
                                                        "-",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3-h),2),"}{",h,"}\\approx",derI,"$$")))
      
    } else {
      output$textoDerivacion1 <-renderText( paste( withMathJax("$$f'(1.3) \\approx \\frac{f(1.3+",h,")-f(1.3)}{",h,"} \\approx$$$$\\frac{ f(",1.3+h,")-f(",1.3,")}{",h,"}\\approx \\frac{ ",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3+h),2),
                                                       "-",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3),2),"}{",h,"}\\approx",derD,"$$")))
      
    }
    
    
    xDer <- seq(1, 1.5, 0.05)
    xDerdf <- as.data.frame(xDer)
    slopeDer <- function(x) {
      xincrementDer <- 0.2
      yincrementDer <-
        ((x+0.1)^(exp(1)-(x+0.1)^2)*exp((-(x+0.1))^2)*(-2*(x+0.1)*ln(x+0.1)+(1/(x+0.1))))-((x-0.1)^(exp(1)-(x-0.1)^2)*exp((-(x-0.1))^2)*(-2*(x-0.1)*ln(x-0.1)+(1/(x-0.1))))
      slopeDer <- yincrementDer / xincrementDer
      return(slopeDer)
    }
    
    interceptDer <- function(x) {
      interceptDer <- (x)^(exp(1)-(x)^2)*exp((-(x))^2)*(-2*(x)*ln(x)+(1/(x))) - slopeDer(x) * x
      yvaluesDer <- slopeDer(x) * xDer + interceptDer
      xDerdf <- as.data.frame(cbind(xDerdf, yvaluesDer))
      return(interceptDer)
    } 
    
    centG <- ggplot(xDerdf, aes(x = xDer, y =  xDer^(exp(1)-xDer^2)*exp((-xDer)^2)*(-2*xDer*ln(xDer)+(1/xDer))))+
      geom_line(aes(x = 1.3,colour = "a"),color="#f2540c", size = 1,linetype = 5) +
      geom_line(aes(x = 1.3+h,colour = "a+h"),color="#f2540c", size = 1,linetype = 5)+
      geom_line(aes(x = 1.3-h,colour = "a-h"),color="#f2540c", size = 1,linetype = 5) +
      geom_abline(
        slope = slopeDer(1.3),
        intercept = interceptDer(1.3),
        color = "black",
        size = 1
      )+annotate(
        "text",
        x = 1.1,
        y = 4,
        label = "f'(a)",
        colour = "black",
        size=6,
        vjust = 1.5
      )+annotate(
        "text",
        x = 1.05,
        y = 2.5,
        label = "f(x)",
        colour = "black",
        size=6,
        vjust = 1.5
      )  +annotate(
        "text",
        x = 1.3,
        y = -6.25,
        label = "a",
        colour = "black",
        size=6,
        vjust = 1.5
      ) + annotate(
        "text",
        x =1.3+h+0.01,
        y =  -6.25,
        label = "a+h",
        colour = "black",
        size=6,
        vjust = 1.5
      )+ annotate(
        "text",
        x =1.3-h-0.01,
        y =  -6.25,
        label = "a-h",
        colour = "black",
        size=6,
        vjust = 1.5
      ) +labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xDerdf) +
      geom_line(color = "blue", size = 1)
    
    izqG <- ggplot(xDerdf, aes(x = xDer, y =  xDer^(exp(1)-xDer^2)*exp((-xDer)^2)*(-2*xDer*ln(xDer)+(1/xDer))))+
      geom_line(aes(x = 1.3,colour = "a"),color="#f2540c", size = 1,linetype = 5) +
      geom_line(aes(x = 1.3-h,colour = "a-h"),color="#f2540c", size = 1,linetype = 5) +
      geom_abline(
        slope = slopeDer(1.3),
        intercept = interceptDer(1.3),
        color = "black",
        size = 1
      )+annotate(
        "text",
        x = 1.1,
        y = 4,
        label = "f'(a)",
        colour = "black",
        size=6,
        vjust = 1.5
      )+annotate(
        "text",
        x = 1.05,
        y = 2.5,
        label = "f(x)",
        colour = "black",
        size=6,
        vjust = 1.5
      )  +annotate(
        "text",
        x = 1.3,
        y = -6.25,
        label = "a",
        colour = "black",
        size=6,
        vjust = 1.5
      ) + annotate(
        "text",
        x =1.3-h-0.01,
        y =  -6.25,
        label = "a-h",
        colour = "black",
        size=6,
        vjust = 1.5
      ) +labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xDerdf) +
      geom_line(color = "blue", size = 1) 
    
    derG <- ggplot(xDerdf, aes(x = xDer, y =  xDer^(exp(1)-xDer^2)*exp((-xDer)^2)*(-2*xDer*ln(xDer)+(1/xDer))))+
      geom_line(aes(x = 1.3,colour = "a"),color="#f2540c", size = 1,linetype = 5) +
      geom_line(aes(x = 1.3+h,colour = "a+h"),color="#f2540c", size = 1,linetype = 5)+
      geom_abline(
        slope = slopeDer(1.3),
        intercept = interceptDer(1.3),
        color = "black",
        size = 1
      )+annotate(
        "text",
        x = 1.1,
        y = 4,
        label = "f '(a)",
        colour = "black",
        size=6,
        vjust = 1.5
      )+annotate(
        "text",
        x = 1.05,
        y = 2.5,
        label = "f(x)",
        colour = "black",
        size=6,
        vjust = 1.5
      )  +annotate(
        "text",
        x = 1.3,
        y = -6.25,
        label = "a",
        colour = "black",
        size=6,
        vjust = 1.5
      ) + annotate(
        "text",
        x =1.3+h+0.01,
        y =  -6.25,
        label = "a+h",
        colour = "black",
        size=6,
        vjust = 1.5
      )+labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xDerdf) +
      geom_line(color = "blue", size = 1)
    
    
    if(input$metodosDerivacion1== "Centrada" ){
      output$graficaDerivacion1 <- renderPlot(centG)
    }else if(input$metodosDerivacion1== "Descentrada izquierda" ){
      output$graficaDerivacion1 <- renderPlot(izqG)
    } else {
      output$graficaDerivacion1 <- renderPlot(derG)
    }
  }
  )
  
  #DERIVACIÓN SENCILLA
  
  observeEvent(ignoreInit = FALSE, c(
    input$metodosDerivacion2, input$h2
  ),{
    req(input$h2)
    
    derC2<- derivC(function(x) 
      (x^2+3*x-2)^4 ,1.3,input$h2)
    derD2<- derivD(function(x)
      (x^2+3*x-2)^4,1.3,input$h2)
    derI2<- derivI(function(x)
      (x^2+3*x-2)^4,1.3,input$h2)
    h2<-input$h2
    
    
    if(input$metodosDerivacion2== "Centrada" ){
      output$textoDerivacion2 <- renderText( paste( withMathJax("$$f'(1.3) \\approx \\frac{f(1.3+",h2,")-f(1.3-",h2,")}{2·",h2,"} \\approx$$ $$\\frac{ f(",1.3+h2,")-f(",1.3-h2,")}{",2*h2,"}\\approx \\frac{ ",round(feval(function(x)(x^2+3*x-2)^4,1.3+h2),2),"-",round(feval(function(x)(x^2+3*x-2)^4,1.3-h2),2),"}{",2*h2,"}\\approx",derC2,"$$")))
      
    }else if(input$metodosDerivacion2== "Descentrada izquierda" ){
      output$textoDerivacion2 <- renderText( paste( withMathJax("$$f'(1.3) \\approx \\frac{f(1.3)-f(1.3-",h2,")}{",h2,"} \\approx$$$$\\frac{ f(",1.3,")-f(",1.3-h2,")}{",h2,"}\\approx \\frac{ ",round(feval(function(x)(x^2+3*x-2)^4,1.3),2),"-",round(feval(function(x)(x^2+3*x-2)^4,1.3-h2),2),"}{",h2,"}\\approx",derI2,"$$")))
      
    } else {
      output$textoDerivacion2 <-renderText( paste( withMathJax("$$f'(1.3) \\approx \\frac{f(1.3+",h2,")-f(1.3)}{",h2,"} \\approx$$$$\\frac{ f(",1.3+h2,")-f(",1.3,")}{",h2,"}\\approx \\frac{ ",round(feval(function(x)(x^2+3*x-2)^4,1.3+h2),2),"-",round(feval(function(x)(x^2+3*x-2)^4,1.3),2),"}{",h2,"}\\approx",derD2,"$$")))
      
    }
    
    output$analiticoDerivacion<-renderText(paste( withMathJax("Resolución de forma analítica: $$f'(x) = (8x+12)(x^2+3x-2)^3$$ $$f'(1.3) = (8·1.3+12)(1.3^2+3x-2)^3=1036.409$$ ")))
    
    xDer2 <- seq(1, 1.5, 0.05)
    xDer2df <- as.data.frame(xDer2)
    slopeDer2 <- function(x) {
      xincrementDer2 <- 0.2
      yincrementDer2 <-
        (((x+0.1)^2+3*(x+0.1)-2)^4)-(((x-0.1)^2+3*(x-0.1)-2)^4)
      slopeDers <- yincrementDer2 / xincrementDer2
      return(slopeDers)
    }
    
    interceptDer2 <- function(x) {
      interceptDers <- ((x^2+3*x-2)^4) - slopeDer2(x) * x
      yvaluesDer2 <- slopeDer2(x) * xDer2 + interceptDers
      xDer2df <- as.data.frame(cbind(xDer2df, yvaluesDer2))
      return(interceptDers)
    } 
    
    centG2 <- ggplot(xDer2df, aes(x = xDer2, y =  ((xDer2^2+3*xDer2-2)^4)))+
      geom_line(aes(x = 1.3,colour = "a"),color="#f2540c", size = 1,linetype = 5) +
      geom_line(aes(x = 1.3+h2,colour = "a+h"),color="#f2540c", size = 1,linetype = 5)+
      geom_line(aes(x = 1.3-h2,colour = "a-h"),color="#f2540c", size = 1,linetype = 5) +
      geom_abline(
        slope = slopeDer2(1.3),
        intercept = interceptDer2(1.3),
        color = "black",
        size = 1
      )+annotate(
        "text",
        x = 1.5,
        y = 350,
        label = "f '(a)",
        colour = "black",
        size=6,
        vjust = 1.5
      )+annotate(
        "text",
        x = 1.5,
        y = 450,
        label = "f(x)",
        colour = "black",
        size=6,
        vjust = 1.5
      )  +annotate(
        "text",
        x = 1.3,
        y = 10,
        label = "a",
        colour = "black",
        size=6,
        vjust = 1.5
      ) + annotate(
        "text",
        x =1.3+h2+0.01,
        y =  10,
        label = "a+h",
        colour = "black",
        size=6,
        vjust = 1.5
      )+ annotate(
        "text",
        x =1.3-h2-0.01,
        y =  10,
        label = "a-h",
        colour = "black",
        size=6,
        vjust = 1.5
      ) +labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xDer2df) +
      geom_line(color = "blue", size = 1)
    
    izqG2 <- ggplot(xDer2df, aes(x = xDer2, y =  ((xDer2^2+3*xDer2-2)^4)))+
      geom_line(aes(x = 1.3,colour = "a"),color="#f2540c", size = 1,linetype = 5) +
      geom_line(aes(x = 1.3-h2,colour = "a-h"),color="#f2540c", size = 1,linetype = 5) +
      geom_abline(
        slope = slopeDer2(1.3),
        intercept = interceptDer2(1.3),
        color = "black",
        size = 1
      )+annotate(
        "text",
        x = 1.5,
        y = 350,
        label = "f'(a)",
        colour = "black",
        size=6,
        vjust = 1.5
      )+annotate(
        "text",
        x = 1.5,
        y = 450,
        label = "f(x)",
        colour = "black",
        size=6,
        vjust = 1.5
      )  +annotate(
        "text",
        x = 1.3,
        y = 10,
        label = "a",
        colour = "black",
        size=6,
        vjust = 1.5
      ) +  annotate(
        "text",
        x =1.3-h2-0.01,
        y =  10,
        label = "a-h",
        colour = "black",
        size=6,
        vjust = 1.5
      ) +labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xDer2df) +
      geom_line(color = "blue", size = 1)
    
    derG2 <- ggplot(xDer2df, aes(x = xDer2, y =  ((xDer2^2+3*xDer2-2)^4)))+
      geom_line(aes(x = 1.3,colour = "a"),color="#f2540c", size = 1,linetype = 5) +
      geom_line(aes(x = 1.3+h2,colour = "a+h"),color="#f2540c", size = 1,linetype = 5)+
      geom_abline(
        slope = slopeDer2(1.3),
        intercept = interceptDer2(1.3),
        color = "black",
        size = 1
      )+annotate(
        "text",
        x = 1.5,
        y = 350,
        label = "f'(a)",
        colour = "black",
        size=6,
        vjust = 1.5
      )+annotate(
        "text",
        x = 1.5,
        y = 450,
        label = "f(x)",
        colour = "black",
        size=6,
        vjust = 1.5
      )  +annotate(
        "text",
        x = 1.3,
        y = 10,
        label = "a",
        colour = "black",
        size=6,
        vjust = 1.5
      ) + annotate(
        "text",
        x =1.3+h2+0.01,
        y =  10,
        label = "a+h",
        colour = "black",
        size=6,
        vjust = 1.5
      ) +labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xDer2df) +
      geom_line(color = "blue", size = 1)
    
    
    if(input$metodosDerivacion2== "Centrada" ){
      output$graficaDerivacion2 <- renderPlot(centG2)
    }else if(input$metodosDerivacion2== "Descentrada izquierda" ){
      output$graficaDerivacion2 <- renderPlot(izqG2)
    } else {
      output$graficaDerivacion2 <- renderPlot(derG2)
    }
    
  }
  )
  
  
  ##INTEGRACIÓN NEWTONCOTES Cerrada ejemplo1----------------------------
  
  observeEvent(input$ReglasIntegracionNewtonC1,{
    
    output$DatosIntegracionNewtonC1<-NULL
    output$TablaIntegracionNewtonC1<-NULL
    output$textIntegracionNewtonC1<-NULL
    output$resultIntegracionNewtonC1<-NULL
    
    output$botonTrapecioNewtonC1<-NULL
    output$botonSimpsonNewtonC1<-NULL
    output$botonSimpson3NewtonC1<-NULL
    output$graficaIntegracionNewtonC1<-NULL
    
    if(input$ReglasIntegracionNewtonC1== "Regla del trapecio"){
      output$botonTrapecioNewtonC1 <- renderUI({
        actionButton("botonTrapecioNewtonC1", label = ">>" ,style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
      })
      
    }else if(input$ReglasIntegracionNewtonC1== "Regla de Simpson" ){
      output$botonSimpsonNewtonC1 <- renderUI({
        actionButton("botonSimpsonNewtonC1", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
      })
      
    } else {
      output$botonSimpson3NewtonC1 <- renderUI({
        actionButton("botonSimpson3NewtonC1", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
      })
    }
    
  })
  observeEvent(input$botonTrapecioNewtonC1,{
    trapecioNCC<- Rtrapecio(function(x) 
      (x-2)/(x-1)^(1/2) ,2, 4)
    xIntegracionNewtonC1 <-seq(1.5, 5, 0.05)
    xIntegracionNewtonC1df <- as.data.frame(xIntegracionNewtonC1)
    
    if(input$botonTrapecioNewtonC1==1){
      output$DatosIntegracionNewtonC1 <- renderText( paste( withMathJax("Sabiendo que \\(x_0 = a = 2, x_n = b = 4\\) y \\(n= 1\\)  $$h =\\frac{b-a}{n}= \\frac{4 -2}{1}=2$$")))
      grafTrapC <- ggplot(xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, y =  (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)))+
        annotate(
          "text",
          x = 5,
          y = 10,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonC1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonC1 <- renderPlot(grafTrapC)
    }
    if(input$botonTrapecioNewtonC1==2){
      output$TablaIntegracionNewtonC1 <- renderText( paste( withMathJax("  $$f(x_0)= \\frac{2-2}{\\sqrt{(2-1)}}=0$$
                                                                        $$f(x_1)= \\frac{4-2}{\\sqrt{(4-1)}}=\\frac{2\\sqrt{3}}{3}$$")))
      
      grafTrapC <- ggplot(xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, y =  (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)))+
        geom_line(aes(x = 2,colour = "a"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "b"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 5,
          y = 10,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 2,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =4,
          y =  0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                    x = 3.07, xend = 3.87, y = -0.25, yend = -0.25,
                    arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 2.93, xend = 2.13, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 3,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonC1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonC1 <- renderPlot(grafTrapC)
    }
    if(input$botonTrapecioNewtonC1==3){
      output$textIntegracionNewtonC1 <- renderText( paste( withMathJax("$$  \\int_{2}^{4}\\frac{x-2}{\\sqrt(x-1)}dx  \\approx \\frac{2}{2}(f(2)+f(4))\\approx$$ ")))
    }
    if(input$botonTrapecioNewtonC1==4){
      output$resultIntegracionNewtonC1 <- renderText( paste( withMathJax("$$1·\\frac{2\\sqrt{3}}{3}\\approx", trapecioNCC,"$$")))
      
      grafTrapC <- ggplot(xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, y =  (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)))+
        geom_ribbon(data = xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, xmax =4, xmin = 2), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, ymax =Inf , ymin =  (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)), fill = "white",alpha=0.7)+
        geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
        
        geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+ 
        geom_line(aes(x = 2,colour = "a"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "b"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 5,
          y = 10,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 2,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =4,
          y =  0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.07, xend = 3.87, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 2.93, xend = 2.13, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 3,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonC1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonC1 <- renderPlot(grafTrapC)
      
    }
    
  })
  
  observeEvent(input$botonSimpsonNewtonC1,{
    simpsonNCC<- Rsimpson(function(x)
      (x-2)/(x-1)^(1/2) ,2, 4)
    xIntegracionNewtonC1 <-  seq(1.5, 5, 0.05)
    xIntegracionNewtonC1df <- as.data.frame(xIntegracionNewtonC1)
    
    if(input$botonSimpsonNewtonC1==1){
      output$DatosIntegracionNewtonC1 <- renderText( paste( withMathJax("Sabiendo que \\( x_i=x_0+ih, x_0 = a = 2, x_n = b = 4\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n}= \\frac{4 -2}{2}=1$$")))
      grafSimpsonC <- ggplot(xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, y =   (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)))+
        annotate(
          "text",
          x = 5,
          y = 10,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC1df) +
        geom_line(color = "blue", size = 1)
      
      output$graficaIntegracionNewtonC1 <- renderPlot(grafSimpsonC)
    }
    if(input$botonSimpsonNewtonC1==2){
      output$TablaIntegracionNewtonC1 <- renderText( paste( withMathJax("$$
                                                                  \\begin{array}{|c|c|c|c|c|}
                                                                              i  & 0 & 1 &2   \\\\
                                                                   \\hline
                                                                         x_i & 2& 3&4 \\\\
                                                                 
                                                                       f(x_i) & 0 & \\frac{\\sqrt{2}}{2}&\\frac{2\\sqrt{3}}{3}  \\\\
                                                                                               
                                                           \\end{array}$$")))
      grafSimpsonC <- ggplot(xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, y =   (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)))+
        geom_line(aes(x = 2,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 3,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x2"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 5,
          y = 10,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 2,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =3,
          y =  0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+ annotate(
          "text",
          x =4,
          y =  0,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 2.6, xend = 2.9, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 2.4, xend = 2.1, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 2.5,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.6, xend = 3.9, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 3.4, xend = 3.1, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 3.5,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC1df) +
        geom_line(color = "blue", size = 1)
      
      output$graficaIntegracionNewtonC1 <- renderPlot(grafSimpsonC)
      
    }
    if(input$botonSimpsonNewtonC1==3){
      output$textIntegracionNewtonC1 <- renderText( paste( withMathJax("$$\\int_{2}^{4}\\frac{x-2}{\\sqrt(x-1)}dx  \\approx \\frac{1}{3}(f(2)+4f(3)+f(4))\\approx$$ ")))
    }
    if(input$botonSimpsonNewtonC1==4){
      output$resultIntegracionNewtonC1 <- renderText( paste( withMathJax("$$\\frac{1}{3}·(0+4·\\frac{\\sqrt{2}}{2}+\\frac{2\\sqrt{3}}{3}) \\approx", simpsonNCC,"$$")))
      grafSimpsonC <- ggplot(xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, y =   (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)))+
        geom_ribbon(data = xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, xmax =4 , xmin = 2), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, ymax =Inf , ymin =  (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)), fill = "white",alpha=0.7)+
        geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
        
        geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
        
        geom_line(aes(x = 2,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 3,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x2"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 5,
          y = 10,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 2,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =3,
          y =  0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+ annotate(
          "text",
          x =4,
          y =  0,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 2.6, xend = 2.9, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 2.4, xend = 2.1, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 2.5,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.6, xend = 3.9, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 3.4, xend = 3.1, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 3.5,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC1df) +
        geom_line(color = "blue", size = 1)
      
      output$graficaIntegracionNewtonC1 <- renderPlot(grafSimpsonC)
      
    }
    
  })
  
  
  observeEvent(input$botonSimpson3NewtonC1,{
    
    simpson3octNCC<- Rsimpson3oct(function(x)
      (x-2)/(x-1)^(1/2) ,2, 4)
    xIntegracionNewtonC1 <- seq(1.5, 5, 0.05)
    xIntegracionNewtonC1df <- as.data.frame(xIntegracionNewtonC1)
    
    if(input$botonSimpson3NewtonC1==1){
      output$DatosIntegracionNewtonC1 <- renderText( paste( withMathJax("Sabiendo que \\( x_i=x_0+ih, x_0 = a = 2, x_n = b = 4\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n}= \\frac{4 -2}{3}=\\frac{2}{3}$$")))
      
      grafSimpson3C <- ggplot(xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, y =   (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)))+
         annotate(
          "text",
          x = 5,
          y = 10,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC1df) +
        geom_line(color = "blue", size = 1)
      
      output$graficaIntegracionNewtonC1 <- renderPlot(grafSimpson3C)
    }
    if(input$botonSimpson3NewtonC1==2){
      output$TablaIntegracionNewtonC1 <- renderText( paste( withMathJax(" $$\\begin{array}{|c|c|c|c|c|}
                                                                        i & 0 & 1 & 2 &3  \\\\
                                                                      \\hline
                                                                      x_i & 2&\\frac{8}{3} &\\frac{10}{3}&4 \\\\
                                                                      
                                                                      f(x_i) &0&  0.52  & 0.87 &\\frac{2\\sqrt{3}}{3}  \\\\
                                                                      
                                                                      \\end{array}$$
                                                                       ")))
      
      grafSimpson3C <- ggplot(xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, y =   (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)))+
        geom_line(aes(x = 2,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = (8)/3,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = (10)/3,colour = "x2"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x3"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 5,
          y = 10,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 2,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =(8)/3,
          y =  0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5 
        )+ annotate(
          "text",
          x =(10)/3,
          y =  0,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+ annotate(
          "text",
          x =4,
          y =  0,
          label = "x3",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 2.45, xend = 2.55, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.15, "cm")))+
        annotate("segment", 
                 x = 2.24, xend = 2.1, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.15, "cm")))+
        annotate(
          "text",
          x = 2.34,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.1, xend = 3.2, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.15, "cm")))+
        annotate("segment", 
                 x = 2.9, xend = 2.8, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.15, "cm")))+
        annotate(
          "text",
          x = 3,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.8, xend = 3.9, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.15, "cm")))+
        annotate("segment", 
                 x = 3.6, xend = 3.5, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.15, "cm")))+
        annotate(
          "text",
          x = 3.7,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=4,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC1df) +
        geom_line(color = "blue", size = 1)
      
      
      output$graficaIntegracionNewtonC1 <- renderPlot(grafSimpson3C)
      
    }
    if(input$botonSimpson3NewtonC1==3){
      output$textIntegracionNewtonC1 <- renderText( paste( withMathJax("$$\\int_{2}^{4}\\frac{x-2}{\\sqrt(x-1)}dx  \\approx\\frac{3·2}{8·3}(f(2)+3f(\\frac{8}{3})+3f(\\frac{10}{3})+f(4))\\approx$$")))
    }
    if(input$botonSimpson3NewtonC1==4){
      output$resultIntegracionNewtonC1 <- renderText( paste( withMathJax("$$\\frac{1}{4}(0+3·0.52+3·0.87+\\frac{2\\sqrt{3}}{3} )\\approx", simpson3octNCC,"$$")))
      
      grafSimpson3C <- ggplot(xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, y =   (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)))+
        geom_ribbon(data = xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, xmax =4 , xmin = 2), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xIntegracionNewtonC1df, aes(x = xIntegracionNewtonC1, ymax =Inf , ymin =  (xIntegracionNewtonC1-2)/(xIntegracionNewtonC1-1)^(1-2)), fill = "white",alpha=0.7)+
        geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
        
        geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
        geom_line(aes(x = 2,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = (8)/3,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = (10)/3,colour = "x2"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x3"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 5,
          y = 10,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 2,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =(8)/3,
          y =  0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5 
        )+ annotate(
          "text",
          x =(10)/3,
          y =  0,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+ annotate(
          "text",
          x =4,
          y =  0,
          label = "x3",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 2.45, xend = 2.55, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.15, "cm")))+
        annotate("segment", 
                 x = 2.24, xend = 2.1, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.15, "cm")))+
        annotate(
          "text",
          x = 2.34,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.1, xend = 3.2, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.15, "cm")))+
        annotate("segment", 
                 x = 2.9, xend = 2.8, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.15, "cm")))+
        annotate(
          "text",
          x = 3,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.8, xend = 3.9, y = -0.25, yend = -0.25,
                   arrow=arrow(length=unit(0.15, "cm")))+
        annotate("segment", 
                 x = 3.6, xend = 3.5, y = -0.25, yend = -0.25,
                 arrow=arrow(length=unit(0.15, "cm")))+
        annotate(
          "text",
          x = 3.7,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=4,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC1df) +
        geom_line(color = "blue", size = 1)
      
      
      output$graficaIntegracionNewtonC1 <- renderPlot(grafSimpson3C)
      
    }
  })
  
  ##INTEGRACIÓN NEWTONCOTES cerrada ejemplo2------------------------------------
  
  observeEvent(input$ReglasIntegracionNewtonC2,{
    
    output$DatosIntegracionNewtonC2<-NULL
    output$TablaIntegracionNewtonC2<-NULL
    output$textIntegracionNewtonC2<-NULL
    output$resultIntegracionNewtonC2<-NULL
    output$AnaliticoIntegracionNewtonC2<-NULL
    
    output$botonTrapecioNewtonC2<-NULL
    output$botonSimpsonNewtonC2<-NULL
    output$botonSimpson3NewtonC2<-NULL
    output$graficaIntegracionNewtonC2<-NULL
    
    if(input$ReglasIntegracionNewtonC2== "Regla del trapecio"){
      output$botonTrapecioNewtonC2 <- renderUI({
        actionButton("botonTrapecioNewtonC2", label = ">>" ,style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
      })
      
    }else if(input$ReglasIntegracionNewtonC2== "Regla de Simpson" ){
      output$botonSimpsonNewtonC2 <- renderUI({
        actionButton("botonSimpsonNewtonC2", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
      })
      
    } else {
      output$botonSimpson3NewtonC2 <- renderUI({
        actionButton("botonSimpson3NewtonC2", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
      })
    }
    
  })
  observeEvent(input$botonTrapecioNewtonC2,{
    trapecio<- Rtrapecio(function(x) 
      x*(x^2+9)^(1/2) ,0, 4 )
    xIntegracionNewtonC2 <- seq(-1, 5, 0.05)
    xIntegracionNewtonC2df <- as.data.frame(xIntegracionNewtonC2)
    
    if(input$botonTrapecioNewtonC2==1){
      output$DatosIntegracionNewtonC2 <- renderText( paste( withMathJax("Sabiendo que \\(x_0 = a = 0, x_n = b = 4\\) y \\(n= 1\\)  $$h =\\frac{b-a}{n}= \\frac{4-0}{1}=4$$")))
      grafTrapS <- ggplot(xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, y =  (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))))+
           annotate(
          "text",
          x = 5,
          y = 25,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC2df) +
        geom_line(color = "blue", size = 1) 
      
      
      output$graficaIntegracionNewtonC2 <- renderPlot(grafTrapS)
    }
    if(input$botonTrapecioNewtonC2==2){
      output$TablaIntegracionNewtonC2 <- renderText( paste( withMathJax("  $$f(x_0)= 0·{\\sqrt{(0^2+9)}}=0$$
                                                                        $$f(x_1)= 4·{\\sqrt{(4^2+9)}}=20$$")))
      
      grafTrapS <- ggplot(xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, y =  (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))))+
       geom_line(aes(x = 0,colour = "a"),color="black", size = 1,linetype = 5) +
        geom_line(aes(x = 4,colour = "b"),color="black", size = 1,linetype = 5)+
        annotate(
          "text",
          x = 5,
          y = 25,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 0,
          y = -5,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =4,
          y =  -5,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                  x = 2.2, xend = 3.8, y = -1, yend = -1,
                  arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 1.8, xend = 0.2, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 2,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC2df) +
        geom_line(color = "blue", size = 1) 
      
      
      output$graficaIntegracionNewtonC2 <- renderPlot(grafTrapS)
    }
    if(input$botonTrapecioNewtonC2==3){
      output$textIntegracionNewtonC2 <- renderText( paste( withMathJax("$$\\int_{0}^{4} x\\sqrt{(x^2+9)} dx \\approx\\frac{4}{2}(f(0)+f(4))\\approx$$ ")))
    }
    if(input$botonTrapecioNewtonC2==4){
      output$resultIntegracionNewtonC2 <- renderText( paste( withMathJax("$$\\frac{1}{2}·(0+20)\\approx", trapecio,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
      
      grafTrapS <- ggplot(xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, y =  (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))))+
        geom_ribbon(data = xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, xmax =4 , xmin = 0), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, ymax =Inf , ymin = (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))), fill = "white",alpha=0.7)+ 
        geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
        
        geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
        geom_line(aes(x = 0,colour = "a"),color="black", size = 1,linetype = 5) +
        geom_line(aes(x = 4,colour = "b"),color="black", size = 1,linetype = 5)+
        annotate(
          "text",
          x = 5,
          y = 25,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 0,
          y = -5,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =4,
          y =  -5,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 2.2, xend = 3.8, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 1.8, xend = 0.2, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 2,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC2df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonC2 <- renderPlot(grafTrapS)
      
    }
    if(input$botonTrapecioNewtonC2==5){
      output$AnaliticoIntegracionNewtonC2<-renderText(paste( withMathJax("$$\\int_{0}^{4} x\\sqrt{(x^2+9)} dx = \\frac{(x^2+9)^{\\frac{3}{2}}}{3} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
                                     $$\\int_{0}^{4} f(x)dx = F(4)-F(0) =\\frac{(4^2+9)^{\\frac{3}{2}}}{3} - \\frac{(0^2+9)^{\\frac{3}{2}}}{3} = 32.67 $$ ")))
      
    }
  })
  
  observeEvent(input$botonSimpsonNewtonC2,{
    simpson<- Rsimpson(function(x)
      x*(x^2+9)^(1/2),0, 4 )
    xIntegracionNewtonC2 <-   seq(-1, 5, 0.05)
    xIntegracionNewtonC2df <- as.data.frame(xIntegracionNewtonC2)
    
    if(input$botonSimpsonNewtonC2==1){
      output$DatosIntegracionNewtonC2 <- renderText( paste( withMathJax("Sabiendo que \\( x_i=x_0+ih, x_0 = a = 0, x_n = b = 4\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n}= \\frac{4-0}{2}=2$$")))
      grafSimpsonS <- ggplot(xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, y =  (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))))+
          annotate(
          "text",
          x = 5,
          y = 25,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC2df) +
        geom_line(color = "blue", size = 1)
      
      output$graficaIntegracionNewtonC2 <- renderPlot(grafSimpsonS)
    }
    if(input$botonSimpsonNewtonC2==2){
      output$TablaIntegracionNewtonC2 <- renderText( paste( withMathJax("$$
                                                                  \\begin{array}{|c|c|c|c|c|}
                                                                              i  & 0 & 1 &2   \\\\
                                                                   \\hline
                                                                         x_i & 0& 2&4 \\\\
                                                                 
                                                                       f(x_i) & 0 & 2\\sqrt{13}&20 \\\\
                                                                                               
                                                           \\end{array}$$")))
      grafSimpsonS <- ggplot(xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, y =  (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))))+
        geom_line(aes(x = 0,colour = "x0"),color="black", size = 1,linetype = 5) +
        geom_line(aes(x = 2,colour = "x1"),color="black", size = 1,linetype = 5) +
        geom_line(aes(x = 4,colour = "x2"),color="black", size = 1,linetype = 5)+
        annotate(
          "text",
          x = 5,
          y = 25,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 0,
          y = -5,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =2,
          y =  -5,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+ annotate(
          "text",
          x =4,
          y =  -5,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                  x = 1.2, xend = 1.8, y = -1, yend = -1,
                  arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 0.8, xend = 0.2, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 1,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.2, xend = 3.8, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 2.8, xend = 2.2, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 3,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC2df) +
        geom_line(color = "blue", size = 1)
      
      output$graficaIntegracionNewtonC2 <- renderPlot(grafSimpsonS)
      
    }
    if(input$botonSimpsonNewtonC2==3){
      output$textIntegracionNewtonC2 <- renderText( paste( withMathJax("  $$\\int_{0}^{4}x\\sqrt{(x^2+9)} dx \\approx \\frac{2}{3}(f(0)+4f(2)+f(4))\\approx$$ ")))
    }
    if(input$botonSimpsonNewtonC2==4){
      output$resultIntegracionNewtonC2 <- renderText( paste( withMathJax("$$\\frac{2}{3}·(0+4·(2\\sqrt{13})+20) \\approx", simpson,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
      grafSimpsonS <- ggplot(xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, y =  (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))))+
        geom_ribbon(data = xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, xmax =4 , xmin = 0), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, ymax =Inf , ymin = (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))), fill = "white",alpha=0.7) +
        geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
        
        geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
        geom_line(aes(x = 0,colour = "x0"),color="black", size = 1,linetype = 5) +
        geom_line(aes(x = 2,colour = "x1"),color="black", size = 1,linetype = 5) +
        geom_line(aes(x = 4,colour = "x2"),color="black", size = 1,linetype = 5)+
        annotate(
          "text",
          x = 5,
          y = 25,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 0,
          y = -5,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =2,
          y =  -5,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+ annotate(
          "text",
          x =4,
          y =  -5,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 1.2, xend = 1.8, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 0.8, xend = 0.2, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 1,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.2, xend = 3.8, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 2.8, xend = 2.2, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 3,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC2df) +
        geom_line(color = "blue", size = 1)
      
      output$graficaIntegracionNewtonC2 <- renderPlot(grafSimpsonS)
      
    }
    if(input$botonSimpsonNewtonC2==5){
      output$AnaliticoIntegracionNewtonC2<-renderText(paste( withMathJax("$$\\int_{0}^{4} x\\sqrt{(x^2+9)} dx = \\frac{(x^2+9)^{\\frac{3}{2}}}{3} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
                                     $$\\int_{0}^{4} f(x)dx = F(4)-F(0) =\\frac{(4^2+9)^{\\frac{3}{2}}}{3} - \\frac{(0^2+9)^{\\frac{3}{2}}}{3} = 32.67 $$ ")))
      
    }
    
  })
  
  
  observeEvent(input$botonSimpson3NewtonC2,{
    
    simpson3oct<- Rsimpson3oct(function(x)
      x*(x^2+9)^(1/2) ,0, 4 )
    xIntegracionNewtonC2 <-  seq(-1, 5, 0.05)
    xIntegracionNewtonC2df <- as.data.frame(xIntegracionNewtonC2)
    
    if(input$botonSimpson3NewtonC2==1){
      output$DatosIntegracionNewtonC2 <- renderText( paste( withMathJax("Sabiendo que \\( x_i=x_0+ih, x_0 = a = 0, x_n = b = 4\\) y \\(n= 3\\)  $$h =\\frac{b-a}{n}= \\frac{4-0}{3}=\\frac{4}{3}$$")))
      
      
      grafSimpson3S <- ggplot(xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, y =  (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))))+
         annotate(
          "text",
          x = 5,
          y = 25,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC2df) +
        geom_line(color = "blue", size = 1)  
      
      
      output$graficaIntegracionNewtonC2 <- renderPlot(grafSimpson3S)
    }
    if(input$botonSimpson3NewtonC2==2){
      output$TablaIntegracionNewtonC2 <- renderText( paste( withMathJax(" $$\\begin{array}{|c|c|c|c|c|}
                                                                        i & 0 & 1 & 2 &3  \\\\
                                                                      \\hline
                                                                      x_i & 0&\\frac{4}{3} &\\frac{8}{3}&4 \\\\
                                                                      
                                                                      f(x_i) &0&  4.38  & 10.7 &20  \\\\
                                                                      
                                                                      \\end{array}$$
                                                                       ")))
      
      
      grafSimpson3S <- ggplot(xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, y =  (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))))+
        geom_line(aes(x = 0,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 1.34,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 2.7,colour = "x2"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x3"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 5,
          y = 25,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 0,
          y = -5,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =1.34,
          y =  -5,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+ annotate(
          "text",
          x =2.7,
          y =  -5,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate(
          "text",
          x =4,
          y =  -5,
          label = "x3",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 0.77, xend = 1.2, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 0.57, xend = 0.1, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 0.67,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+annotate("segment", 
                   x = 2.1, xend = 2.5, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 1.9, xend = 1.5, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 2,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.5, xend = 3.85, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 3.2, xend = 2.8, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 10/3,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC2df) +
        geom_line(color = "blue", size = 1)  
      
      
      
      output$graficaIntegracionNewtonC2 <- renderPlot(grafSimpson3S)
      
    }
    if(input$botonSimpson3NewtonC2==3){
      output$textIntegracionNewtonC2 <- renderText( paste( withMathJax(" $$\\int_{0}^{4} x\\sqrt{(x^2+9)} dx \\approx \\frac{3·4}{8·3}(f(0)+3f(\\frac{4}{3})+3f(\\frac{8}{3})+f(4))\\approx$$")))
    }
    if(input$botonSimpson3NewtonC2==4){
      output$resultIntegracionNewtonC2 <- renderText( paste( withMathJax("$$\\frac{1}{2}(0+3·4.38+3·10.7+20 )\\approx", simpson3oct,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
      
      
      grafSimpson3S <- ggplot(xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, y =  (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))))+
        geom_ribbon(data = xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, xmax =4 , xmin = 0), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xIntegracionNewtonC2df, aes(x = xIntegracionNewtonC2, ymax =Inf , ymin = (xIntegracionNewtonC2*(xIntegracionNewtonC2^2+9)^(1/2))), fill = "white",alpha=0.7)+
        geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
        
        geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
        geom_line(aes(x = 0,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 1.34,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 2.7,colour = "x2"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x3"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 5,
          y = 25,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 0,
          y = -5,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =1.34,
          y =  -5,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+ annotate(
          "text",
          x =2.7,
          y =  -5,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate(
          "text",
          x =4,
          y =  -5,
          label = "x3",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 0.77, xend = 1.2, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 0.57, xend = 0.1, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 0.67,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+annotate("segment", 
                   x = 2.1, xend = 2.5, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 1.9, xend = 1.5, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 2,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.5, xend = 3.85, y = -1, yend = -1,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 3.2, xend = 2.8, y = -1, yend = -1,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 10/3,
          y = -0.0025,
          label = "h",
          colour = "black",
          size=5,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +theme_bw()+
        geom_line(data = xIntegracionNewtonC2df) +
        geom_line(color = "blue", size = 1)  
      
      output$graficaIntegracionNewtonC2 <- renderPlot(grafSimpson3S)
      
    }
    if(input$botonSimpson3NewtonC2==5){
      output$AnaliticoIntegracionNewtonC2<-renderText(paste( withMathJax("$$\\int_{0}^{4} x\\sqrt{(x^2+9)} dx = \\frac{(x^2+9)^{\\frac{3}{2}}}{3} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
                                     $$\\int_{0}^{4} f(x)dx = F(4)-F(0) =\\frac{(4^2+9)^{\\frac{3}{2}}}{3} - \\frac{(0^2+9)^{\\frac{3}{2}}}{3} = 32.67 $$ ")))
      
    }
  })
  

    ##---------------------INTEGRACIÓN NEWTONCOTES Abierta ejemplo1-----------
  observeEvent(input$ReglasIntegracionNewtonA1,{
    
    output$DatosIntegracionNewtonA1<-NULL
    output$TablaIntegracionNewtonA1<-NULL
    output$textIntegracionNewtonA1<-NULL
    output$resultIntegracionNewtonA1<-NULL
    output$AnaliticoIntegracionNewtonA1<-NULL
    
    output$botonPMNewtonA1<-NULL
    output$boton2PNewtonA1<-NULL
    output$boton3PNewtonA1<-NULL
    output$graficaIntegracionNewtonA1<-NULL
    
    if(input$ReglasIntegracionNewtonA1== "Regla del punto medio"){
      output$botonPMNewtonA1 <- renderUI({
        actionButton("botonPMNewtonA1", label = ">>" ,style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
      })
      
    }else if(input$ReglasIntegracionNewtonA1== "Fórmula de 2 puntos" ){
      output$boton2PNewtonA1 <- renderUI({
        actionButton("boton2PNewtonA1", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
      })
      
    } else {
      output$boton3PNewtonA1 <- renderUI({
        actionButton("boton3PNewtonA1", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
      })
    }
    
  })
  observeEvent(input$botonPMNewtonA1,{
    PuntoMedioNA1<- RPuntoMedio(function(x) 
      (x^3)/((x^2+1)*(x^2+4)*(x^2+9)) ,3, 5)
    xIntegracionNewtonA1 <- seq(0, 6, 0.05)
    xIntegracionNewtonA1df <- as.data.frame(xIntegracionNewtonA1)
    
    if(input$botonPMNewtonA1==1){
      output$DatosIntegracionNewtonA1 <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a=3, x_{n+1} = b=5\\) y \\(n= 0\\)  $$h =\\frac{b-a}{n+2}= \\frac{5 -3}{2}=1$$")))
      grafPuntoM1 <- ggplot(xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, y = (xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9))))+
        labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonA1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonA1 <- renderPlot(grafPuntoM1)
    }
    if(input$botonPMNewtonA1==2){
      output$TablaIntegracionNewtonA1 <- renderText( paste( withMathJax(" $$x_0= a+h = 4$$ $$f(x_0)= \\frac{4^3}{(4^2+1)(4^2+4)(4^2+9)}=0.008$$")))
      grafPuntoM1 <- ggplot(xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, y = (xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9))))+
       geom_line(aes(x = 3,colour = "a"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 5,colour = "b"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 1,
          y = 0.014,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 3,
          y = -0,
          label = "x(-1)",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +annotate(
          "text",
          x = 4,
          y =  -0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =5,
          y =   -0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.57, xend = 3.87, y = 0.002, yend = 0.002,
                   arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 3.43, xend = 3.13, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 3.5,
          y = 0.0025,
          label = "h",
          colour = "black",
          size=6,
          vjust = 1.5
        )+
        annotate("segment", 
                 x = 4.57, xend = 4.87, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 4.43, xend = 4.13, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 4.5,
          y = 0.0025,
          label ="h",
          colour = "black",
          size=6,
          vjust = 1.5
        )  +
        labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonA1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonA1 <- renderPlot(grafPuntoM1)
    }
    if(input$botonPMNewtonA1==3){
      output$textIntegracionNewtonA1 <- renderText( paste( withMathJax("  $$\\int_{3}^{5}\\frac{x^3}{(x^2+1)(x^2+4)(x^2+9)}dx \\approx 2·1(f(4))\\approx $$ ")))
    }
    if(input$botonPMNewtonA1==4){
      output$resultIntegracionNewtonA1 <- renderText( paste( withMathJax("$$2·0.008\\approx", PuntoMedioNA1,"$$")))
      grafPuntoM1 <- ggplot(xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, y = (xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9))))+
        geom_ribbon(data = xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1 , xmax =5 , xmin = 3), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, ymax =Inf , ymin = ((xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9)))), fill = "white",alpha=0.7)+
        geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
        
        geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
        geom_line(aes(x = 3,colour = "a"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 5,colour = "b"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 1,
          y = 0.014,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 3,
          y = -0,
          label = "x(-1)",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +annotate(
          "text",
          x = 4,
          y =  -0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =5,
          y =   -0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                  x = 3.57, xend = 3.87, y = 0.002, yend = 0.002,
                  arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 3.43, xend = 3.13, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 3.5,
          y = 0.0025,
          label = "h",
          colour = "black",
          size=6,
          vjust = 1.5
        )+
        annotate("segment", 
                 x = 4.57, xend = 4.87, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate("segment", 
                 x = 4.43, xend = 4.13, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.3, "cm")))+
        annotate(
          "text",
          x = 4.5,
          y = 0.0025,
          label ="h",
          colour = "black",
          size=6,
          vjust = 1.5
        )  +
        labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonA1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonA1 <- renderPlot(grafPuntoM1)
      
    }
    
  })
  
  observeEvent(input$boton2PNewtonA1,{
    abierta2pNA1<- fabierta2P(function(x)
      (x^3)/((x^2+1)*(x^2+4)*(x^2+9)) ,3, 5)
    xIntegracionNewtonA1 <- seq(0, 6, 0.05)
    xIntegracionNewtonA1df <- as.data.frame(xIntegracionNewtonA1)
    
    if(input$boton2PNewtonA1==1){
      output$DatosIntegracionNewtonA1 <- renderText( paste( withMathJax("Sabiendo que \\(x_0= a+h, x_1= b-h,x_{-1} = a=3, x_{n+1} = b=5\\) y \\(n= 1\\)  $$h =\\frac{b-a}{n+2}= \\frac{5 -3}{1+2}=\\frac{2}{3}$$ ")))
      graf2P1<-ggplot(xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, y = (xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9))))+
        labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonA1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonA1 <- renderPlot(graf2P1)
    }
    if(input$boton2PNewtonA1==2){
      output$TablaIntegracionNewtonA1 <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
       i  & 0 & 1    \\\\
       \\hline
       x_i & \\frac{11}{3} &\\frac{13}{3} \\\\
     
      f(x_i) & 0.0087 & 0.0065  \\\\
                                   
      \\end{array}$$")))
      
      graf2P1<-ggplot(xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, y = (xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9))))+
        geom_line(aes(x = 3,colour = "a"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 11/3,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 13/3,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 5,colour = "b"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 1,
          y = 0.014,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 3,
          y = 0,
          label = "x(-1)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 11/3,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +annotate(
          "text",
          x = 13/3,
          y = 0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =5,
          y =  0,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.37, xend = 3.62, y = 0.002, yend = 0.002,
                   arrow=arrow(length=unit(0.25, "cm")))+
        annotate("segment", 
                 x = 3.23, xend = 3.02, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate(
          "text",
          x = 3.3,
          y = 0.0025,
          label = "h",
          colour = "black",
          size=6,
          vjust = 1.5
        )+
        annotate("segment", 
                 x = 4.07, xend = 4.28, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate("segment", 
                 x = 3.93, xend = 3.68, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate(
          "text",
          x = 4,
          y = 0.0025,
          label ="h",
          colour = "black",
          size=6,
          vjust = 1.5
        )  +
        annotate("segment", 
                 x = 4.7, xend = 4.92, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate("segment", 
                 x = 4.55, xend = 4.35, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate(
          "text",
          x = 4.666,
          y = 0.0025,
          label ="h",
          colour = "black",
          size=6,
          vjust = 1.5
        ) +labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonA1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonA1 <- renderPlot(graf2P1)
      
    }
    if(input$boton2PNewtonA1==3){
      output$textIntegracionNewtonA1 <- renderText( paste( withMathJax("$$\\int_{3}^{5}\\frac{x^3}{(x^2+1)(x^2+4)(x^2+9)}dx  \\approx \\frac{3\\frac{2}{3}}{2}(f(\\frac{11}{3})+f(\\frac{13}{3}))\\approx$$ ")))
    }
    if(input$boton2PNewtonA1==4){
      output$resultIntegracionNewtonA1 <- renderText( paste( withMathJax("$$1·(0.0087+0.0065)\\approx", abierta2pNA1,"$$")))
      graf2P1<-ggplot(xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, y = (xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9))))+
        geom_ribbon(data = xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1 , xmax =5 , xmin = 3), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, ymax =Inf , ymin = ((xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9)))), fill = "white",alpha=0.7)+
        geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
        
        geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
        geom_line(aes(x = 3,colour = "a"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 11/3,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 13/3,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 5,colour = "b"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 1,
          y = 0.014,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 3,
          y = 0,
          label = "x(-1)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 11/3,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +annotate(
          "text",
          x = 13/3,
          y = 0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =5,
          y =  0,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        )+annotate("segment", 
                   x = 3.37, xend = 3.62, y = 0.002, yend = 0.002,
                   arrow=arrow(length=unit(0.25, "cm")))+
        annotate("segment", 
                 x = 3.23, xend = 3.02, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate(
          "text",
          x = 3.3,
          y = 0.0025,
          label = "h",
          colour = "black",
          size=6,
          vjust = 1.5
        )+
        annotate("segment", 
                 x = 4.07, xend = 4.28, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate("segment", 
                 x = 3.93, xend = 3.68, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate(
          "text",
          x = 4,
          y = 0.0025,
          label ="h",
          colour = "black",
          size=6,
          vjust = 1.5
        ) +
        annotate("segment", 
                 x = 4.7, xend = 4.92, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate("segment", 
                 x = 4.55, xend = 4.35, y = 0.002, yend = 0.002,
                 arrow=arrow(length=unit(0.25, "cm")))+
        annotate(
          "text",
          x = 4.666,
          y = 0.0025,
          label ="h",
          colour = "black",
          size=6,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonA1df) +
        geom_line(color = "blue", size = 1) 
      
       
      
      output$graficaIntegracionNewtonA1 <- renderPlot(graf2P1)
      
    }
    
  })
  
  
  observeEvent(input$boton3PNewtonA1,{
    abierta3pNA1<- fabierta3P(function(x)
      (x^3)/((x^2+1)*(x^2+4)*(x^2+9)) ,3, 5)
    xIntegracionNewtonA1 <- seq(0, 6, 0.05)
    xIntegracionNewtonA1df <- as.data.frame(xIntegracionNewtonA1)
    
    if(input$boton3PNewtonA1==1){
      output$DatosIntegracionNewtonA1 <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a, x_{n+1}=b,x_0= a+h,x_i= x_0+ih, x_n= b-h\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n+2}= \\frac{5 -3}{2+2}=\\frac{1}{2}$$")))
      graf3P1<-ggplot(xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, y = (xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9))))+
      labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonA1df) +
        geom_line(color = "blue", size = 1) 
      output$graficaIntegracionNewtonA1 <- renderPlot(graf3P1)
    }
    if(input$boton3PNewtonA1==2){
      output$TablaIntegracionNewtonA1 <- renderText( paste( withMathJax(" $$\\begin{array}{|c|c|c|c|c|}
                                                                        i &-1 & 0 & 1 & 2 &3  \\\\
                                                                      \\hline
                                                                      x_i & 3&\\frac{7}{2} &4 &\\frac{9}{2}&5 \\\\
                                                                      
                                                                      f(x_i) &0.0115&  0.0094  & 0.0075 & 0.006 &0.0049  \\\\
                                                                      
                                                                      \\end{array}$$
                                                                        
                                                                        En este caso, al igual que en los anteriores ejemplos, el tamaño de los subintervalos es h. ")))
      
      
      graf3P1<-ggplot(xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, y = (xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9))))+
        geom_line(aes(x = 3,colour = "a"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 7/2,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 9/2,colour = "x2"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 5,colour = "b"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 1,
          y = 0.014,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 3,
          y = 0,
          label = "x(-1)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 7/2,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +annotate(
          "text",
          x = 4,
          y = 0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +annotate(
          "text",
          x = 9/2,
          y = 0,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =5,
          y =  0,
          label = "x3",
          colour = "black",
          size=4,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonA1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonA1 <- renderPlot(graf3P1)
      
    }
    if(input$boton3PNewtonA1==3){
      output$textIntegracionNewtonA1 <- renderText( paste( withMathJax("$$\\int_{3}^{5}\\frac{x^3}{(x^2+1)(x^2+4)(x^2+9)}dx  \\approx \\frac{4\\frac{1}{2}}{3}(2f(\\frac{7}{2})-f(4)+2f(\\frac{9}{2}))\\approx$$")))
    }
    if(input$boton3PNewtonA1==4){
      output$resultIntegracionNewtonA1 <- renderText( paste( withMathJax("$$\\frac{2}{3}(2·0.0094-0.0075+2·0.006)\\approx", abierta3pNA1,"$$")))
      
      
      graf3P1<-ggplot(xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, y = (xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9))))+
        geom_ribbon(data = xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1 , xmax =5 , xmin = 3), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xIntegracionNewtonA1df, aes(x = xIntegracionNewtonA1, ymax =Inf , ymin = ((xIntegracionNewtonA1^3)/((xIntegracionNewtonA1^2+1)*(xIntegracionNewtonA1^2+4)*(xIntegracionNewtonA1^2+9)))), fill = "white",alpha=0.7)+
        geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
        
        geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
        geom_line(aes(x = 3,colour = "a"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 7/2,colour = "x0"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 4,colour = "x1"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 9/2,colour = "x2"),color="black", size = 1,linetype = 3) +
        geom_line(aes(x = 5,colour = "b"),color="black", size = 1,linetype = 3)+
        annotate(
          "text",
          x = 1,
          y = 0.014,
          label = "f(x)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 3,
          y = 0,
          label = "x(-1)",
          colour = "black",
          size=4,
          vjust = 1.5
        )  +annotate(
          "text",
          x = 7/2,
          y = 0,
          label = "x0",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +annotate(
          "text",
          x = 4,
          y = 0,
          label = "x1",
          colour = "black",
          size=4,
          vjust = 1.5
        ) +annotate(
          "text",
          x = 9/2,
          y = 0,
          label = "x2",
          colour = "black",
          size=4,
          vjust = 1.5
        ) + annotate(
          "text",
          x =5,
          y =  0,
          label = "x3",
          colour = "black",
          size=4,
          vjust = 1.5
        )+labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xIntegracionNewtonA1df) +
        geom_line(color = "blue", size = 1) 
      
      output$graficaIntegracionNewtonA1 <- renderPlot(graf3P1)
      
    }
 
  })
  
  
    
    ##---------------------INTEGRACIÓN NEWTONCOTES Abierta ejemplo2------------
    observeEvent(input$ReglasIntegracionNewtonA2,{
      
      output$DatosIntegracionNewtonA2<-NULL
      output$TablaIntegracionNewtonA2<-NULL
      output$textIntegracionNewtonA2<-NULL
      output$resultIntegracionNewtonA2<-NULL
      output$AnaliticoIntegracionNewtonA2<-NULL
      
      output$botonPMNewtonA2<-NULL
      output$boton2PNewtonA2<-NULL
      output$boton3PNewtonA2<-NULL
      output$graficaIntegracionNewtonA2<-NULL
      
      if(input$ReglasIntegracionNewtonA2== "Regla del punto medio"){
        output$botonPMNewtonA2 <- renderUI({
          actionButton("botonPMNewtonA2", label = ">>" ,style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
        
      }else if(input$ReglasIntegracionNewtonA2== "Fórmula de 2 puntos" ){
        output$boton2PNewtonA2 <- renderUI({
          actionButton("boton2PNewtonA2", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
        
      } else {
        output$boton3PNewtonA2 <- renderUI({
          actionButton("boton3PNewtonA2", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
      }
      
    })
    observeEvent(input$botonPMNewtonA2,{
      PuntoMedioNAS<- RPuntoMedio(function(x) 
        cos(x)*sin(x) ,0, pi/2)
      xIntegracionNewtonA2 <- seq(-0.25, 1.75, 0.05)
      xIntegracionNewtonA2df <- as.data.frame(xIntegracionNewtonA2)
      
      if(input$botonPMNewtonA2==1){
        output$DatosIntegracionNewtonA2 <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a=0, x_{n+1} = b=\\frac{\\pi}{2}\\) y \\(n= 0\\)  $$h =\\frac{b-a}{n+2}= \\frac{\\frac{\\pi}{2}-0}{2}=\\frac{\\pi}{4}$$")))
        grafPuntoMS<-ggplot(xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, y = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)))+
          annotate(
            "text",
            x = -0.2,
            y = -0.125,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )+labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionNewtonA2df) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionNewtonA2 <- renderPlot(grafPuntoMS)
         }
      if(input$botonPMNewtonA2==2){
        output$TablaIntegracionNewtonA2 <- renderText( paste( withMathJax("$$x_0= a+h =\\frac{\\pi}{4}$$ $$f(x_0)= cos(\\frac{\\pi}{4})sin(\\frac{\\pi}{4})=\\frac{1}{2}$$")))
        
        grafPuntoMS <- ggplot(xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, y = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)))+
          geom_line(aes(x = 0,colour = "a"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/4,colour = "x0"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/2,colour = "b"),color="black", size = 1,linetype = 3)+
          annotate(
            "text",
            x = -0.2,
            y = -0.125,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0,
            y = -0.25,
            label = "x(-1)",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = pi/4,
            y = -0.25,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =pi/2,
            y =  -0.25,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +
          annotate("segment", 
                   x = 0.47, xend = 0.713, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.25, "cm")))+
          annotate("segment", 
                   x = 0.32, xend = 0.073, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.25, "cm")))+
          annotate(
            "text",
            x = 0.393,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          )+
          annotate("segment", 
                   x = 1.25, xend = 1.5, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.25, "cm")))+
          annotate("segment", 
                   x = 1.1, xend = 0.85, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.25, "cm")))+
          annotate(
            "text",
            x = 1.179,
            y = -0.002,
            label ="h",
            colour = "black",
            size=6,
            vjust = 1.5
          )  + labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionNewtonA2df) +
          geom_line(color = "blue", size = 1) 
        
        
        output$graficaIntegracionNewtonA2 <- renderPlot(grafPuntoMS)
      }
      if(input$botonPMNewtonA2==3){
        output$textIntegracionNewtonA2 <- renderText( paste( withMathJax("  $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx  \\approx 2·\\frac{\\pi}{4}(f(\\frac{\\pi}{4}))\\approx $$ ")))
      }
      if(input$botonPMNewtonA2==4){
        output$resultIntegracionNewtonA2 <- renderText( paste( withMathJax("$$\\frac{\\pi}{2}·\\frac{1}{2} \\approx", PuntoMedioNAS,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
        grafPuntoMS <- ggplot(xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, y = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)))+
          geom_ribbon(data = xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, xmax =pi/2 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, ymax =Inf , ymin = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)), fill = "white",alpha=0.7)+
          geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
          geom_line(aes(x = 0,colour = "a"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/4,colour = "x0"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/2,colour = "b"),color="black", size = 1,linetype = 3)+
          geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
          
          annotate(
            "text",
            x = -0.2,
            y = -0.125,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0,
            y = -0.25,
            label = "x(-1)",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = pi/4,
            y = -0.25,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =pi/2,
            y =  -0.25,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +
          annotate("segment", 
                   x = 0.47, xend = 0.713, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.25, "cm")))+
          annotate("segment", 
                   x = 0.32, xend = 0.073, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.25, "cm")))+
          annotate(
            "text",
            x = 0.393,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          )+
          annotate("segment", 
                   x = 1.25, xend = 1.5, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.25, "cm")))+
          annotate("segment", 
                   x = 1.1, xend = 0.85, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.25, "cm")))+
          annotate(
            "text",
            x = 1.179,
            y = -0.002,
            label ="h",
            colour = "black",
            size=6,
            vjust = 1.5
          )  + labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionNewtonA2df) +
          geom_line(color = "blue", size = 1) 
        
        
        output$graficaIntegracionNewtonA2 <- renderPlot(grafPuntoMS)
        
      }
      if(input$botonPMNewtonA2==5){
      output$AnaliticoIntegracionNewtonA2<-renderText(paste( withMathJax(" $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx = \\frac{sen(x)^2}{2} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
    $$\\int_{0}^{\\frac{\\pi}{2}} f(x)dx = F(\\frac{\\pi}{2})-F(0) =  \\frac{sen(\\frac{\\pi}{2})^2}{2} - \\frac{sen(0)^2}{2}= 0.5 $$ ")))
      }
    })
    
    observeEvent(input$boton2PNewtonA2,{
      abierta2pNAS<- fabierta2P(function(x)
        cos(x)*sin(x) ,0, pi/2)
      xIntegracionNewtonA2 <- seq(-0.25, 1.75, 0.05)
      xIntegracionNewtonA2df <- as.data.frame(xIntegracionNewtonA2)
      
      if(input$boton2PNewtonA2==1){
        output$DatosIntegracionNewtonA2 <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a, x_{n+1} = b , x_0= a+h ,x_n= b-h \\) y \\(n= 1\\)  $$h =\\frac{b-a}{n+2}= \\frac{\\frac{\\pi}{2}-0}{1+2}=\\frac{\\pi}{6}$$ ")))
        graf2PS<-ggplot(xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, y = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)))+
          annotate(
            "text",
            x = -0.2,
            y = -0.125,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )+labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionNewtonA2df) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionNewtonA2 <- renderPlot(graf2PS)
        }
      if(input$boton2PNewtonA2==2){
        output$TablaIntegracionNewtonA2 <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
       i & -1 & 0 & 1 & 2   \\\\
       \\hline
       x_i & 0 & \\frac{\\pi}{6} &\\frac{\\pi}{3} & \\frac{\\pi}{2}\\\\
     
      f(x_i) & 0  & \\frac{\\sqrt{3}}{4} & \\frac{\\sqrt{3}}{4} & 0 \\\\
                                   
      \\end{array}$$")))
        
        graf2PS<-ggplot(xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, y =cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)))+
         geom_line(aes(x = 0,colour = "a"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/6,colour = "x0"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/3,colour = "x1"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/2,colour = "b"),color="black", size = 1,linetype = 3)+
          annotate(
            "text",
            x = -0.2,
            y = -0.125,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0,
            y = -0.25,
            label = "x(-1)",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = pi/6,
            y = -0.25,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =pi/3,
            y =  -0.25,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =pi/2,
            y =  -0.25,
            label = "x2",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +
          annotate("segment", 
                   x = 0.332, xend = 0.482, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate("segment", 
                   x = 0.192, xend = 0.042, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate(
            "text",
            x = 0.262,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          )+
          annotate("segment", 
                   x = 0.855, xend = 1, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate("segment", 
                   x = 0.72, xend = 0.565, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate(
            "text",
            x = 0.785,
            y = -0.002,
            label ="h",
            colour = "black",
            size=6,
            vjust = 1.5
          )  +
          annotate("segment", 
                   x = 1.38, xend = 1.53, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate("segment", 
                   x = 1.24, xend = 1.09, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate(
            "text",
            x = 1.31,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          ) +labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionNewtonA2df) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionNewtonA2 <- renderPlot(graf2PS)
        
      }
      if(input$boton2PNewtonA2==3){
        output$textIntegracionNewtonA2 <- renderText( paste( withMathJax("$$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx \\approx \\frac{3\\frac{\\pi}{6}}{2}(f(\\frac{\\pi}{6})+f(\\frac{\\pi}{3}))\\approx $$ ")))
      }
      if(input$boton2PNewtonA2==4){
        output$resultIntegracionNewtonA2 <- renderText( paste( withMathJax("$$\\frac{\\pi}{4}(\\frac{\\sqrt{3}}{4}+\\frac{\\sqrt{3}}{4})\\approx", abierta2pNAS,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
        graf2PS<-ggplot(xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, y =cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)))+
          geom_ribbon(data = xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, xmax =pi/2 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, ymax =Inf , ymin = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)), fill = "white",alpha=0.7)+
          geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
          
          geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
          geom_line(aes(x = 0,colour = "a"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/6,colour = "x0"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/3,colour = "x1"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/2,colour = "b"),color="black", size = 1,linetype = 3)+
          annotate(
            "text",
            x = -0.2,
            y = -0.125,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0,
            y = -0.25,
            label = "x(-1)",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = pi/6,
            y = -0.25,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =pi/3,
            y =  -0.25,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =pi/2,
            y =  -0.25,
            label = "x2",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +
          annotate("segment", 
                   x = 0.332, xend = 0.482, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate("segment", 
                   x = 0.192, xend = 0.042, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate(
            "text",
            x = 0.262,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          )+
          annotate("segment", 
                   x = 0.855, xend = 1, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate("segment", 
                   x = 0.72, xend = 0.565, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate(
            "text",
            x = 0.785,
            y = -0.002,
            label ="h",
            colour = "black",
            size=6,
            vjust = 1.5
          )  +
          annotate("segment", 
                   x = 1.38, xend = 1.53, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate("segment", 
                   x = 1.24, xend = 1.09, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.15, "cm")))+
          annotate(
            "text",
            x = 1.31,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          ) +labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionNewtonA2df) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionNewtonA2 <- renderPlot(graf2PS)
         
      }
      if(input$boton2PNewtonA2==5){
        output$AnaliticoIntegracionNewtonA2<-renderText(paste( withMathJax(" $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx = \\frac{sen(x)^2}{2} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
    $$\\int_{0}^{\\frac{\\pi}{2}} f(x)dx = F(\\frac{\\pi}{2})-F(0) =  \\frac{sen(\\frac{\\pi}{2})^2}{2} - \\frac{sen(0)^2}{2}= 0.5 $$ ")))
      }
    })
    
    
    observeEvent(input$boton3PNewtonA2,{
      abierta3pNAS<- fabierta3P(function(x)
        cos(x)*sin(x) ,0, pi/2)
      xIntegracionNewtonA2 <- seq(-0.25, 1.75, 0.05)
      xIntegracionNewtonA2df <- as.data.frame(xIntegracionNewtonA2)
      
      if(input$boton3PNewtonA2==1){
        output$DatosIntegracionNewtonA2 <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a, x_{n+1} = b,x_0= a+h,x_i= x_0+ih, x_n= b-h\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n+2}= \\frac{\\frac{\\pi}{2}-0}{2+2}=\\frac{\\pi}{8}$$")))
        graf3PS<-ggplot(xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, y = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)))+
            annotate(
            "text",
            x = -0.2,
            y = -0.125,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )+labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionNewtonA2df) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionNewtonA2 <- renderPlot(graf3PS)
      }
      if(input$boton3PNewtonA2==2){
        output$TablaIntegracionNewtonA2 <- renderText( paste( withMathJax(" $$\\begin{array}{|c|c|c|c|c|}
                                                                        i & -1 & 0 & 1 & 2 & 3  \\\\
                                                                      \\hline
                                                                      x_i & 0 & \\frac{\\pi}{8} &\\frac{\\pi}{4} &\\frac{3\\pi}{8} & \\frac{\\pi}{2}\\\\
                                                                      
                                                                      f(x_i) & 0  & 0.354  & \\frac{\\pi}{2} & 0.354 &0 \\\\
                                                                      
                                                                      \\end{array}$$")))
        
        
        graf3PS<-ggplot(xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, y = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)))+
         geom_line(aes(x = 0,colour = "a"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/8,colour = "x0"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/4,colour = "x1"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = (3*pi)/8,colour = "x2"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/2,colour = "b"),color="black", size = 1,linetype = 3)+
          annotate(
            "text",
            x = -0.2,
            y = -0.125,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0,
            y = -0.25,
            label = "x(-1)",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = pi/8,
            y = -0.25,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =pi/4,
            y =  -0.25,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =(3*pi)/8,
            y =  -0.25,
            label = "x2",
            colour = "black",
            size=4,
            vjust = 1.5
          )+ annotate(
            "text",
            x =pi/2,
            y =  -0.25,
            label = "x3",
            colour = "black",
            size=4,
            vjust = 1.5
          )+
          annotate("segment", 
                   x = 0.266, xend = 0.35, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.10, "cm")))+
          annotate("segment", 
                   x = 0.126, xend = 0.026, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.10, "cm")))+
          annotate(
            "text",
            x = 0.196,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          )+
          annotate("segment", 
                   x = 0.658, xend = 0.758, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate("segment", 
                   x = 0.518, xend = 0.418, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate(
            "text",
            x = 0.588,
            y = -0.002,
            label ="h",
            colour = "black",
            size=6,
            vjust = 1.5
          )  +
          annotate("segment", 
                   x = 1.05, xend = 1.15, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate("segment", 
                   x = 0.91, xend = 0.81, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate(
            "text",
            x = 0.98,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          )  +
          annotate("segment", 
                   x = 1.44, xend = 1.54, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate("segment", 
                   x = 1.3, xend = 1.2, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate(
            "text",
            x = 1.37,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          ) +labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionNewtonA2df) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionNewtonA2 <- renderPlot(graf3PS)
        
      }
      if(input$boton3PNewtonA2==3){
        output$textIntegracionNewtonA2 <- renderText( paste( withMathJax("$$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx \\approx\\frac{4\\frac{\\pi}{8}}{3}(2f(\\frac{\\pi}{8})-f(\\frac{\\pi}{4})+2f(\\frac{3\\pi}{8}))\\approx$$")))
      }
      if(input$boton3PNewtonA2==4){
        output$resultIntegracionNewtonA2 <- renderText( paste( withMathJax("$$\\frac{\\pi}{6}(2·0.354-\\frac{1}{2}+2·0.354)\\approx", abierta3pNAS,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
        
        graf3PS<-ggplot(xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, y = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)))+
          geom_ribbon(data = xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, xmax =pi/2 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xIntegracionNewtonA2df, aes(x = xIntegracionNewtonA2, ymax =Inf , ymin = cos(xIntegracionNewtonA2)*sin(xIntegracionNewtonA2)), fill = "white",alpha=0.7)+
          geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
          
          geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
          geom_line(aes(x = 0,colour = "a"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/8,colour = "x0"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/4,colour = "x1"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = (3*pi)/8,colour = "x2"),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = pi/2,colour = "b"),color="black", size = 1,linetype = 3)+
          annotate(
            "text",
            x = -0.2,
            y = -0.125,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +
          annotate(
            "text",
            x = 0,
            y = -0.25,
            label = "x(-1)",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = pi/8,
            y = -0.25,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =pi/4,
            y =  -0.25,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =(3*pi)/8,
            y =  -0.25,
            label = "x2",
            colour = "black",
            size=4,
            vjust = 1.5
          )+ annotate(
            "text",
            x =pi/2,
            y =  -0.25,
            label = "x3",
            colour = "black",
            size=4,
            vjust = 1.5
          )+
          annotate("segment", 
                   x = 0.266, xend = 0.35, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.10, "cm")))+
          annotate("segment", 
                   x = 0.126, xend = 0.026, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.10, "cm")))+
          annotate(
            "text",
            x = 0.196,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          )+
          annotate("segment", 
                   x = 0.658, xend = 0.758, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate("segment", 
                   x = 0.518, xend = 0.418, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate(
            "text",
            x = 0.588,
            y = -0.002,
            label ="h",
            colour = "black",
            size=6,
            vjust = 1.5
          )  +
          annotate("segment", 
                   x = 1.05, xend = 1.15, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate("segment", 
                   x = 0.91, xend = 0.81, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate(
            "text",
            x = 0.98,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          )  +
          annotate("segment", 
                   x = 1.44, xend = 1.54, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate("segment", 
                   x = 1.3, xend = 1.2, y = -0.025, yend = -0.025,
                   arrow=arrow(length=unit(0.1, "cm")))+
          annotate(
            "text",
            x = 1.37,
            y = -0.002,
            label = "h",
            colour = "black",
            size=6,
            vjust = 1.5
          ) +labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionNewtonA2df) +
          geom_line(color = "blue", size = 1) 
        output$graficaIntegracionNewtonA2 <- renderPlot(graf3PS)
       
      }
      if(input$boton3PNewtonA2==5){
        output$AnaliticoIntegracionNewtonA2<-renderText(paste( withMathJax("$$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx = \\frac{sen(x)^2}{2} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
    $$\\int_{0}^{\\frac{\\pi}{2}} f(x)dx = F(\\frac{\\pi}{2})-F(0) =  \\frac{sen(\\frac{\\pi}{2})^2}{2} - \\frac{sen(0)^2}{2}= 0.5 $$ ")))
      }
    })
    
    
    
    
    
    ##---------------------------INTEGRACIÓN COMPUESTA---------------------
    
    
##---------------------Integración compuesta ejemplo 1----------------------
    observeEvent(input$ReglasIntegracionComp1,{
    
      output$DatosIntegracionComp1<-NULL
      output$TablaIntegracionComp1<-NULL
      output$textoIntegracionComp1<-NULL
      output$resultIntegracionComp1<-NULL
      
      output$botonTrapecioComp1<-NULL
      output$botonSimpsonComp1<-NULL
      output$botonPMComp1<-NULL
      output$graficaIntegracionComp1<-NULL
      
      if(input$ReglasIntegracionComp1== "Regla del trapecio" ){
        output$botonTrapecioComp1 <- renderUI({
          actionButton("botonTrapecioComp1", label = ">>" ,style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
        
      }else if(input$ReglasIntegracionComp1== "Regla de Simpson" ){
        output$botonSimpsonComp1 <- renderUI({
          actionButton("botonSimpsonComp1", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
        
      } else {
        output$botonPMComp1 <- renderUI({
          actionButton("botonPMComp1", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
      }
      
    })
    
    observeEvent(input$botonTrapecioComp1,{
      trapecioComp<- RtrapecioC(function(x) 
        abs(x-2)^3*(1-sin(pi*x)),0, 2,4)
      
      xIntegracionComp <- seq(-0.002, 2.5, 0.05)
      xIntegracionCompdf <- as.data.frame(xIntegracionComp)
      
      if(input$botonTrapecioComp1==1){
        output$DatosIntegracionComp1 <- renderText( paste( withMathJax("Sabiendo que \\(x_0=a=0, x_N=b=2\\) y \\(N= 4 \\)   $$h =\\frac{b-a}{N}= \\frac{2 -0}{4}=\\frac{1}{2}$$")))
        grafComp <- ggplot(xIntegracionCompdf, aes(x = xIntegracionComp, y =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))))+
          annotate(
            "text",
            x = 2.4,
            y = 1,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +
          labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionComp1 <- renderPlot(grafComp)
      }
      if(input$botonTrapecioComp1==2){
        output$TablaIntegracionComp1 <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
       i & 0 & 1 &2 & 3 &4 \\\\
       \\hline
       x_i & 0 & 0.5 &1 & 1.5 & 2\\\\
     
      f(x_i) & 8  & 0  & 1 & \\frac{1}{4} &0 \\\\
                                   
      \\end{array}$$")))
        
        grafComp <- ggplot(xIntegracionCompdf, aes(x = xIntegracionComp, y =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))))+
          geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = 0.5),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = 1),color="black", size = 1,linetype = 3)+
          geom_line(aes(x = 1.5),color="black", size = 1,linetype = 3)+
          geom_line(aes(x = 2),color="black", size = 1,linetype = 3)+
          annotate(
            "text",
            x = 2.4,
            y = 1,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0,
            y = 0,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0.5,
            y = 0,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = 1,
            y = 0,
            label = "x2",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = 1.5,
            y = 0,
            label = "x3",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =2,
            y =  0,
            label = "x4",
            colour = "black",
            size=4,
            vjust = 1.5
          )+
          labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionComp1 <- renderPlot(grafComp)
        
      }
      if(input$botonTrapecioComp1==3){
        output$textoIntegracionComp1 <- renderText( paste( withMathJax("$$\\int_{0}^{2}|x-2|^3(1-sen(\\pi x))dx \\approx \\frac{\\frac{1}{2}}{2}(f(0)+2(f(0.5)+f(1)+f(1.5))+f(2))\\approx$$ ")))
        
      }
      if(input$botonTrapecioComp1==4){
        output$resultIntegracionComp1 <- renderText( paste( withMathJax("$$\\frac{1}{4}(8+2(0+1+\\frac{1}{4})+0)\\approx", trapecioComp,"$$")))
        grafComp <- ggplot(xIntegracionCompdf, aes(x = xIntegracionComp, y =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))))+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =0.5 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =1 , xmin = 0.5), fill = "#ff3d70",alpha=0.8)+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =1.5 , xmin = 1), fill = "#d172da",alpha=0.8)+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =2 , xmin = 1.5), fill = "#51a1fa",alpha=0.8)+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp, ymax =Inf , ymin =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))), fill = "white",alpha=0.7)+
          geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
          
          geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
          geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = 0.5),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = 1),color="black", size = 1,linetype = 3)+
          geom_line(aes(x = 1.5),color="black", size = 1,linetype = 3)+
          geom_line(aes(x = 2),color="black", size = 1,linetype = 3)+
          annotate(
            "text",
            x = 2.4,
            y = 1,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0,
            y = 0,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0.5,
            y = 0,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = 1,
            y = 0,
            label = "x2",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = 1.5,
            y = 0,
            label = "x3",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =2,
            y =  0,
            label = "x4",
            colour = "black",
            size=4,
            vjust = 1.5
          )+
          labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionComp1 <- renderPlot(grafComp)
        
      }
      })
      
      observeEvent(input$botonSimpsonComp1,{
      simpsonComp<- RsimpsonC(function(x)
        abs(x-2)^3*(1-sin(pi*x)),0, 2,4)
      
      xIntegracionComp <- seq(-0.002, 2.5, 0.05)
      xIntegracionCompdf <- as.data.frame(xIntegracionComp)
      if(input$botonSimpsonComp1==1){
        output$DatosIntegracionComp1 <- renderText( paste( withMathJax("Sabiendo que \\(x_0=a=0, x_N=b=2\\) y \\(N= 4\\)  $$h =\\frac{b-a}{N}= \\frac{2 -0}{4}=\\frac{1}{2}$$")))
        grafComp <- ggplot(xIntegracionCompdf, aes(x = xIntegracionComp, y =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))))+
          annotate(
            "text",
            x = 2.4,
            y = 1,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +
          labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionComp1 <- renderPlot(grafComp)
      }
      if(input$botonSimpsonComp1==2){
        output$TablaIntegracionComp1 <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
                                                                 i & 0 & 1 &2 & 3 &4 \\\\
                                                                 \\hline
                                                                  x_i & 0 & \\frac{1}{2} &1 & \\frac{3}{2} & 2\\\\
                                                               
                                                                f(x_i) & 8  & 0  & 1 & \\frac{1}{4} &0 \\\\
                                                                                             
                                                                \\end{array}$$")))
        grafComp <- ggplot(xIntegracionCompdf, aes(x = xIntegracionComp, y =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))))+
          geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = 0.5),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = 1),color="black", size = 1,linetype = 3)+
          geom_line(aes(x = 1.5),color="black", size = 1,linetype = 3)+
          geom_line(aes(x = 2),color="black", size = 1,linetype = 3)+
          annotate(
            "text",
            x = 2.4,
            y = 1,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0,
            y = 0,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0.5,
            y = 0,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = 1,
            y = 0,
            label = "x2",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = 1.5,
            y = 0,
            label = "x3",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =2,
            y =  0,
            label = "x4",
            colour = "black",
            size=4,
            vjust = 1.5
          )+
          labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionComp1 <- renderPlot(grafComp)
        
      }
      if(input$botonSimpsonComp1==3){
        output$textoIntegracionComp1 <- renderText( paste( withMathJax("$$\\int_{0}^{2}|x-2|^3(1-sen(\\pi x))dx \\approx\\frac{\\frac{1}{2}}{3}(f(0)+2(f(1))+4(f(\\frac{1}{2})+f(\\frac{3}{2}))+f(2))\\approx $$")))
      }
      
      if(input$botonSimpsonComp1==4){
        output$resultIntegracionComp1 <- renderText( paste( withMathJax("$$\\frac{1}{6}(8+2(1)+4(0+\\frac{1}{4})+0)\\approx ",simpsonComp,"$$")))
        grafComp <- ggplot(xIntegracionCompdf, aes(x = xIntegracionComp, y =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))))+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =0.5 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =1 , xmin = 0.5), fill = "#ff3d70",alpha=0.8)+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =1.5 , xmin = 1), fill = "#d172da",alpha=0.8)+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =2 , xmin = 1.5), fill = "#51a1fa",alpha=0.8)+
          geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp, ymax =Inf , ymin =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))), fill = "white",alpha=0.7)+
          geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
          
          geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
          geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = 0.5),color="black", size = 1,linetype = 3) +
          geom_line(aes(x = 1),color="black", size = 1,linetype = 3)+
          geom_line(aes(x = 1.5),color="black", size = 1,linetype = 3)+
          geom_line(aes(x = 2),color="black", size = 1,linetype = 3)+
          annotate(
            "text",
            x = 2.4,
            y = 1,
            label = "f(x)",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0,
            y = 0,
            label = "x0",
            colour = "black",
            size=4,
            vjust = 1.5
          )  +annotate(
            "text",
            x = 0.5,
            y = 0,
            label = "x1",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = 1,
            y = 0,
            label = "x2",
            colour = "black",
            size=4,
            vjust = 1.5
          ) +annotate(
            "text",
            x = 1.5,
            y = 0,
            label = "x3",
            colour = "black",
            size=4,
            vjust = 1.5
          ) + annotate(
            "text",
            x =2,
            y =  0,
            label = "x4",
            colour = "black",
            size=4,
            vjust = 1.5
          )+
          labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xIntegracionCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaIntegracionComp1 <- renderPlot(grafComp)
      }
      })
      
      observeEvent(input$botonPMComp1,{
        puntoMedioComp<- RPuntoMedioC(function(x)
          abs(x-2)^3*(1-sin(pi*x)),0, 2,4)
        
        xIntegracionComp <- seq(-0.002, 2.5, 0.05)
        xIntegracionCompdf <- as.data.frame(xIntegracionComp)
        if(input$botonPMComp1==1){
          output$DatosIntegracionComp1 <- renderText( paste( withMathJax("Sabiendo que \\(x_0=a=0, x_N=b=2\\) y \\(N= 4\\)  $$h =\\frac{b-a}{N+2}= \\frac{2 -0}{4+2}=\\frac{1}{3}$$")))
          grafComp <- ggplot(xIntegracionCompdf, aes(x = xIntegracionComp, y =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))))+
           annotate(
              "text",
              x = 2.4,
              y = 1,
              label = "f(x)",
              colour = "black",
              size=4,
              vjust = 1.5
            )  +
            labs(x = "X", y = "F(X)") +
            theme_bw()+
            geom_line(data = xIntegracionCompdf) +
            geom_line(color = "blue", size = 1) 
          
          output$graficaIntegracionComp1 <- renderPlot(grafComp)
          }
        if(input$botonPMComp1==2){
          output$TablaIntegracionComp1 <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
                                                                 i & 0 & 1 &2 & 3 &4 \\\\
                                                                 \\hline
                                                                 x_i & 0 & \\frac{1}{2} &1 & \\frac{3}{2} & 2\\\\
                                                               
                                                                f(x_i) & 8  & 0  & 1 & \\frac{1}{4} &0 \\\\
                                                                                             
                                                                \\end{array}$$")))
          grafComp <- ggplot(xIntegracionCompdf, aes(x = xIntegracionComp, y =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))))+
            geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
            geom_line(aes(x = 0.5),color="black", size = 1,linetype = 3) +
            geom_line(aes(x = 1),color="black", size = 1,linetype = 3)+
            geom_line(aes(x = 1.5),color="black", size = 1,linetype = 3)+
            geom_line(aes(x = 2),color="black", size = 1,linetype = 3)+
            annotate(
              "text",
              x = 2.4,
              y = 1,
              label = "f(x)",
              colour = "black",
              size=4,
              vjust = 1.5
            )  +annotate(
              "text",
              x = 0,
              y = 0,
              label = "x0",
              colour = "black",
              size=4,
              vjust = 1.5
            )  +annotate(
              "text",
              x = 0.5,
              y = 0,
              label = "x1",
              colour = "black",
              size=4,
              vjust = 1.5
            ) +annotate(
              "text",
              x = 1,
              y = 0,
              label = "x2",
              colour = "black",
              size=4,
              vjust = 1.5
            ) +annotate(
              "text",
              x = 1.5,
              y = 0,
              label = "x3",
              colour = "black",
              size=4,
              vjust = 1.5
            ) + annotate(
              "text",
              x =2,
              y =  0,
              label = "x4",
              colour = "black",
              size=4,
              vjust = 1.5
            )+
            labs(x = "X", y = "F(X)") +
            theme_bw()+
            geom_line(data = xIntegracionCompdf) +
            geom_line(color = "blue", size = 1) 
          
          output$graficaIntegracionComp1 <- renderPlot(grafComp)
          
        }
        if(input$botonPMComp1==3){
          output$textoIntegracionComp1 <- renderText( paste( withMathJax("$$\\int_{0}^{2}|x-2|^3(1-sen(\\pi x))dx \\approx 2·\\frac{1}{3}(f(0)+f(1)+f(2))\\approx $$")))
        }
        
        if(input$botonPMComp1==4){
          output$resultIntegracionComp1 <- renderText( paste( withMathJax("$$\\frac{2}{3}(8+1+0)\\approx ",puntoMedioComp,"$$")))
          grafComp <- ggplot(xIntegracionCompdf, aes(x = xIntegracionComp, y =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))))+
            geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =0.5 , xmin = 0), fill = "#f2540c",alpha=0.8)+
            geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =1 , xmin = 0.5), fill = "#ff3d70",alpha=0.8)+
            geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =1.5 , xmin = 1), fill = "#d172da",alpha=0.8)+
            geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp , xmax =2 , xmin = 1.5), fill = "#51a1fa",alpha=0.8)+
            geom_ribbon(data = xIntegracionCompdf, aes(x = xIntegracionComp, ymax =Inf , ymin =  abs(xIntegracionComp-2)^3*(1-sin(pi*xIntegracionComp))), fill = "white",alpha=0.7)+
            geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
            
            geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
            geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
            geom_line(aes(x = 0.5),color="black", size = 1,linetype = 3) +
            geom_line(aes(x = 1),color="black", size = 1,linetype = 3)+
            geom_line(aes(x = 1.5),color="black", size = 1,linetype = 3)+
            geom_line(aes(x = 2),color="black", size = 1,linetype = 3)+
            annotate(
              "text",
              x = 2.4,
              y = 1,
              label = "f(x)",
              colour = "black",
              size=4,
              vjust = 1.5
            )  +annotate(
              "text",
              x = 0,
              y = 0,
              label = "x0",
              colour = "black",
              size=4,
              vjust = 1.5
            )  +annotate(
              "text",
              x = 0.5,
              y = 0,
              label = "x1",
              colour = "black",
              size=4,
              vjust = 1.5
            ) +annotate(
              "text",
              x = 1,
              y = 0,
              label = "x2",
              colour = "black",
              size=4,
              vjust = 1.5
            ) +annotate(
              "text",
              x = 1.5,
              y = 0,
              label = "x3",
              colour = "black",
              size=4,
              vjust = 1.5
            ) + annotate(
              "text",
              x =2,
              y =  0,
              label = "x4",
              colour = "black",
              size=4,
              vjust = 1.5
            )+
            labs(x = "X", y = "F(X)") +
            theme_bw()+
            geom_line(data = xIntegracionCompdf) +
            geom_line(color = "blue", size = 1) 
          
          output$graficaIntegracionComp1 <- renderPlot(grafComp)
          
        }
      })
        
##---------------------Integración compuesta ejemplo 2----------------------
        
        observeEvent(input$ReglasIntegracionComp2,{
          
          output$DatosIntegracionComp2<-NULL
          output$TablaIntegracionComp2<-NULL
          output$textoIntegracionComp2<-NULL
          output$resultIntegracionComp2<-NULL
          output$AnaliticoIntegracionComp2<-NULL
          
          output$botonTrapecioComp2<-NULL
          output$botonSimpsonComp2<-NULL
          output$botonPMComp2<-NULL
          output$graficaIntegracionComp2<-NULL
          
          if(input$ReglasIntegracionComp2== "Regla del trapecio" ){
            output$botonTrapecioComp2 <- renderUI({
              actionButton("botonTrapecioComp2", label = ">>" ,style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
            })
            
          }else if(input$ReglasIntegracionComp2== "Regla de Simpson" ){
            output$botonSimpsonComp2 <- renderUI({
              actionButton("botonSimpsonComp2", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
            })
            
          } else {
            output$botonPMComp2 <- renderUI({
              actionButton("botonPMComp2", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
            })
          }
          
        })
        
        observeEvent(input$botonTrapecioComp2,{
          trapecioComp<- RtrapecioC(function(x) 
            (x^3+6*x^2),-6, 0,6)
          
          xIntegracionComp2 <- seq(-6.5, 0.5, 0.05)
          xIntegracionComp2df <- as.data.frame(xIntegracionComp2)
          
          if(input$botonTrapecioComp2==1){
            output$DatosIntegracionComp2 <- renderText( paste( withMathJax("Sabiendo que \\(x_0=a=-6, x_N=b=0\\) y \\(N= 6 \\)   $$h =\\frac{b-a}{N}= \\frac{0 -(-6)}{6}=1$$")))
            grafComp2 <- ggplot(xIntegracionComp2df, aes(x = xIntegracionComp2, y =   (xIntegracionComp2^3+6*xIntegracionComp2^2)))+
              annotate(
                "text",
                x =-6.5,
                y =  -10,
                label = "f(x)",
                colour = "black",
                size=4,
                vjust = 1.5
              )+
              labs(x = "X", y = "F(X)") +
              theme_bw()+
              geom_line(data = xIntegracionComp2df) +
              geom_line(color = "blue", size = 1) 
            
            output$graficaIntegracionComp2 <- renderPlot(grafComp2)
          }
          if(input$botonTrapecioComp2==2){
            output$TablaIntegracionComp2 <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
                                                                 i & 0 & 1 &2 & 3 &4& 5 &6 \\\\
                                                                 \\hline
                                                                 x_i & -6 & -5 &-4 & -3 & -2& -1 &0\\\\
                                                               
                                                                f(x_i) & 0  & 25  & 32 &27 & 16 &5 &0 \\\\
                                                                                             
                                                                \\end{array}$$")))
            
            grafComp2 <- ggplot(xIntegracionComp2df, aes(x = xIntegracionComp2, y =  (xIntegracionComp2^3+6*xIntegracionComp2^2)))+
              geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -1),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -2),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -3),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -4),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -5),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -6),color="black", size = 1,linetype = 3)+
              annotate(
                "text",
                x = -6.5,
                y = -10,
                label = "f(x)",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -6,
                y = -21.5,
                label = "x0",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -5,
                y = -21.5,
                label = "x1",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -4,
                y = -21.5,
                label = "x2",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -3,
                y = -21.5,
                label = "x3",
                colour = "black",
                size=4,
                vjust = 1.5
              ) + annotate(
                "text",
                x =-2,
                y =  -21.5,
                label = "x4",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =-1,
                y =  -21.5,
                label = "x5",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =0,
                y =  -21.5,
                label = "x6",
                colour = "black",
                size=4,
                vjust = 1.5
              )+
              labs(x = "X", y = "F(X)") +
              theme_bw()+
              geom_line(data = xIntegracionComp2df) +
              geom_line(color = "blue", size = 1) 
            
            output$graficaIntegracionComp2 <- renderPlot(grafComp2)
          }
          if(input$botonTrapecioComp2==3){
            output$textoIntegracionComp2 <- renderText( paste( withMathJax("$$\\int_{-6}^{0}x^3+6x^2dx \\approx 
                                                                          $$$$\\frac{1}{2}(f(-6)+2(f(-5)+f(-4)+f(-3)+f(-2)+f(-1))+f(0))\\approx$$ ")))
            
          }
          if(input$botonTrapecioComp2==4){
            output$resultIntegracionComp2 <- renderText( paste( withMathJax("$$\\frac{1}{2}(0+2·(25+32+27+16+5)+0)\\approx", trapecioComp,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
            grafComp2 <- ggplot(xIntegracionComp2df, aes(x = xIntegracionComp2, y =  (xIntegracionComp2^3+6*xIntegracionComp2^2)))+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-5 , xmin = -6), fill = "#f2540c",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-4 , xmin = -5), fill = "#ff4345",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-3 , xmin = -4), fill = "#ff4899",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-2 , xmin = -3), fill = "#d172da",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-1 , xmin = -2), fill = "#8194f9",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =0 , xmin = -1), fill = "#0caaf2",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2, ymax =Inf , ymin =  (xIntegracionComp2^3+6*xIntegracionComp2^2)), fill = "white",alpha=0.7)+
              geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
              
              geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
              geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -1),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -2),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -3),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -4),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -5),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -6),color="black", size = 1,linetype = 3)+
              annotate(
                "text",
                x = -6.5,
                y = -10,
                label = "f(x)",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -6,
                y = -21.5,
                label = "x0",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -5,
                y = -21.5,
                label = "x1",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -4,
                y = -21.5,
                label = "x2",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -3,
                y = -21.5,
                label = "x3",
                colour = "black",
                size=4,
                vjust = 1.5
              ) + annotate(
                "text",
                x =-2,
                y =  -21.5,
                label = "x4",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =-1,
                y =  -21.5,
                label = "x5",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =0,
                y =  -21.5,
                label = "x6",
                colour = "black",
                size=4,
                vjust = 1.5
              )+
              labs(x = "X", y = "F(X)") +
              theme_bw()+
              geom_line(data = xIntegracionComp2df) +
              geom_line(color = "blue", size = 1) 
            
            output$graficaIntegracionComp2 <- renderPlot(grafComp2)
            
          }
          if(input$botonTrapecioComp2==5){
            output$AnaliticoIntegracionComp2<-renderText(paste( withMathJax("$$\\int_{-6}^{0}x^3+6x^2dx = \\frac{(x^3·(x+8))}{4} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
                                     $$\\int_{-6}^{0} f(x)dx = F(0)-F(-6) =\\frac{(0^3·(0+8)}{4} - \\frac{(-6)^3·(-6+8)}{4} = 108 $$ ")))
            
          }
        })
        
        observeEvent(input$botonSimpsonComp2,{
          simpsonComp<- RsimpsonC(function(x)
            (x^3+6*x^2),-6, 0,6)
          
          xIntegracionComp2 <- seq(-6.5, 0.5, 0.05)
          xIntegracionComp2df <- as.data.frame(xIntegracionComp2)
          if(input$botonSimpsonComp2==1){
            output$DatosIntegracionComp2 <- renderText( paste( withMathJax("Sabiendo que \\(x_0=a=-6, x_N=b=0\\) y \\(N= 6\\)  $$h =\\frac{b-a}{N}= \\frac{0 -(-6)}{6}=1$$")))
            grafComp2 <- ggplot(xIntegracionComp2df, aes(x = xIntegracionComp2, y =   (xIntegracionComp2^3+6*xIntegracionComp2^2)))+
              annotate(
                "text",
                x =-6.5,
                y =  -10,
                label = "f(x)",
                colour = "black",
                size=4,
                vjust = 1.5
              )+
              labs(x = "X", y = "F(X)") +
              theme_bw()+
              geom_line(data = xIntegracionComp2df) +
              geom_line(color = "blue", size = 1) 
            
            output$graficaIntegracionComp2 <- renderPlot(grafComp2)
          }
          if(input$botonSimpsonComp2==2){
            output$TablaIntegracionComp2 <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
                                                                           i & 0 & 1 &2 & 3 &4& 5 &6 \\\\
                                                                           \\hline
                                                                           x_i & -6 & -5 &-4 & -3 & -2& -1 &0\\\\
                                                                           
                                                                           f(x_i) & 0  & 25  & 32 &27 & 16 &5 &0 \\\\
                                                                           
                                                                           \\end{array}$$")))
            grafComp2 <- ggplot(xIntegracionComp2df, aes(x = xIntegracionComp2, y =  (xIntegracionComp2^3+6*xIntegracionComp2^2)))+
              geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -1),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -2),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -3),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -4),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -5),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -6),color="black", size = 1,linetype = 3)+
              annotate(
                "text",
                x = -6.5,
                y = -10,
                label = "f(x)",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -6,
                y = -21.5,
                label = "x0",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -5,
                y = -21.5,
                label = "x1",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -4,
                y = -21.5,
                label = "x2",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -3,
                y = -21.5,
                label = "x3",
                colour = "black",
                size=4,
                vjust = 1.5
              ) + annotate(
                "text",
                x =-2,
                y =  -21.5,
                label = "x4",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =-1,
                y =  -21.5,
                label = "x5",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =0,
                y =  -21.5,
                label = "x6",
                colour = "black",
                size=4,
                vjust = 1.5
              )+
              labs(x = "X", y = "F(X)") +
              theme_bw()+
              geom_line(data = xIntegracionComp2df) +
              geom_line(color = "blue", size = 1) 
            
            output$graficaIntegracionComp2 <- renderPlot(grafComp2)
            
          }
          if(input$botonSimpsonComp2==3){
            output$textoIntegracionComp2 <- renderText( paste( withMathJax("$$\\int_{-6}^{0}x^3+6x^2dx \\approx$$$$\\frac{1}{3}(f(-6)+2(f(-4)+f(-2))+4(f(-5)+f(-3)+f(-1))+f(0))\\approx $$")))
          }
          
          if(input$botonSimpsonComp2==4){
            output$resultIntegracionComp2 <- renderText( paste( withMathJax("$$\\frac{1}{3}(0+2(32+16)+4(25+27+5)+0)\\approx ",simpsonComp,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
            grafComp2 <- ggplot(xIntegracionComp2df, aes(x = xIntegracionComp2, y =  (xIntegracionComp2^3+6*xIntegracionComp2^2)))+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-5 , xmin = -6), fill = "#f2540c",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-4 , xmin = -5), fill = "#ff4345",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-3 , xmin = -4), fill = "#ff4899",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-2 , xmin = -3), fill = "#d172da",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-1 , xmin = -2), fill = "#8194f9",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =0 , xmin = -1), fill = "#0caaf2",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2, ymax =Inf , ymin =  (xIntegracionComp2^3+6*xIntegracionComp2^2)), fill = "white",alpha=0.7)+
              geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
              
              geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
              geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -1),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -2),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -3),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -4),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -5),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -6),color="black", size = 1,linetype = 3)+
              annotate(
                "text",
                x = -6.5,
                y = -10,
                label = "f(x)",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -6,
                y = -21.5,
                label = "x0",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -5,
                y = -21.5,
                label = "x1",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -4,
                y = -21.5,
                label = "x2",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -3,
                y = -21.5,
                label = "x3",
                colour = "black",
                size=4,
                vjust = 1.5
              ) + annotate(
                "text",
                x =-2,
                y =  -21.5,
                label = "x4",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =-1,
                y =  -21.5,
                label = "x5",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =0,
                y =  -21.5,
                label = "x6",
                colour = "black",
                size=4,
                vjust = 1.5
              )+
              labs(x = "X", y = "F(X)") +
              theme_bw()+
              geom_line(data = xIntegracionComp2df) +
              geom_line(color = "blue", size = 1) 
            
            output$graficaIntegracionComp2 <- renderPlot(grafComp2)
          } 
          if(input$botonSimpsonComp2==5){
            output$AnaliticoIntegracionComp2<-renderText(paste( withMathJax("$$\\int_{-6}^{0}x^3+6x^2dx = \\frac{(x^3·(x+8))}{4} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
                                     $$\\int_{-6}^{0} f(x)dx = F(0)-F(-6) =\\frac{(0^3·(0+8)}{4} - \\frac{(-6)^3·(-6+8)}{4} = 108 $$ ")))
            
          }
        })
        
        observeEvent(input$botonPMComp2,{
          puntoMedioComp<- RPuntoMedioC(function(x)
            (x^3+6*x^2),-6, 0,6)
          
          xIntegracionComp2 <- seq(-6.5, 0.5, 0.05)
          xIntegracionComp2df <- as.data.frame(xIntegracionComp2)
          if(input$botonPMComp2==1){
            output$DatosIntegracionComp2 <- renderText( paste( withMathJax("Sabiendo que \\(x_0=a=-6, x_N=b=0\\) y \\(N= 6\\)  $$h =\\frac{b-a}{N+2}= \\frac{0 -(-6)}{6+2}=\\frac{3}{4}$$")))
            grafComp2 <- ggplot(xIntegracionComp2df, aes(x = xIntegracionComp2, y =   (xIntegracionComp2^3+6*xIntegracionComp2^2)))+
               annotate(
                "text",
                x =-6.5,
                y =  -10,
                label = "f(x)",
                colour = "black",
                size=4,
                vjust = 1.5
              )+
              labs(x = "X", y = "F(X)") +
              theme_bw()+
              geom_line(data = xIntegracionComp2df) +
              geom_line(color = "blue", size = 1) 
            
            output$graficaIntegracionComp2 <- renderPlot(grafComp2)
          }
          if(input$botonPMComp2==2){
            output$TablaIntegracionComp2 <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
                                                                 i & 0 & 1 &2 & 3 &4& 5 &6 \\\\
                                                                 \\hline
                                                                 x_i & -6 & -5 &-4 & -3 & -2& -1 &0\\\\
                                                               
                                                                f(x_i) & 0  & 25  & 32 &27 & 16 &5 &0 \\\\
                                                                                             
                                                                \\end{array}$$")))
            grafComp2 <- ggplot(xIntegracionComp2df, aes(x = xIntegracionComp2, y =  (xIntegracionComp2^3+6*xIntegracionComp2^2)))+
              geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -1),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -2),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -3),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -4),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -5),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -6),color="black", size = 1,linetype = 3)+
              annotate(
                "text",
                x = -6.5,
                y = -10,
                label = "f(x)",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -6,
                y = -21.5,
                label = "x0",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -5,
                y = -21.5,
                label = "x1",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -4,
                y = -21.5,
                label = "x2",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -3,
                y = -21.5,
                label = "x3",
                colour = "black",
                size=4,
                vjust = 1.5
              ) + annotate(
                "text",
                x =-2,
                y =  -21.5,
                label = "x4",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =-1,
                y =  -21.5,
                label = "x5",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =0,
                y =  -21.5,
                label = "x6",
                colour = "black",
                size=4,
                vjust = 1.5
              )+
              labs(x = "X", y = "F(X)") +
              theme_bw()+
              geom_line(data = xIntegracionComp2df) +
              geom_line(color = "blue", size = 1) 
            
            output$graficaIntegracionComp2 <- renderPlot(grafComp2)
            
          }
          if(input$botonPMComp2==3){
            output$textoIntegracionComp2 <- renderText( paste( withMathJax("$$\\int_{-6}^{0}x^3+6x^2dx \\approx 2·\\frac{3}{4}(f(-6)+f(-4)+f(-2)+f(0))\\approx $$")))
          }
          
          if(input$botonPMComp2==4){
            output$resultIntegracionComp2 <- renderText( paste( withMathJax("$$\\frac{3}{2}(0+32+16+0)\\approx ",puntoMedioComp,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
            grafComp2 <- ggplot(xIntegracionComp2df, aes(x = xIntegracionComp2, y =  (xIntegracionComp2^3+6*xIntegracionComp2^2)))+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-5 , xmin = -6), fill = "#f2540c",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-4 , xmin = -5), fill = "#ff4345",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-3 , xmin = -4), fill = "#ff4899",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-2 , xmin = -3), fill = "#d172da",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =-1 , xmin = -2), fill = "#8194f9",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2 , xmax =0 , xmin = -1), fill = "#0caaf2",alpha=0.8)+
              geom_ribbon(data = xIntegracionComp2df, aes(x = xIntegracionComp2, ymax =Inf , ymin =  (xIntegracionComp2^3+6*xIntegracionComp2^2)), fill = "white",alpha=0.7)+
              geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
              
              geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
              geom_line(aes(x = 0),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -1),color="black", size = 1,linetype = 3) +
              geom_line(aes(x = -2),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -3),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -4),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -5),color="black", size = 1,linetype = 3)+
              geom_line(aes(x = -6),color="black", size = 1,linetype = 3)+
              annotate(
                "text",
                x = -6.5,
                y = -10,
                label = "f(x)",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -6,
                y = -21.5,
                label = "x0",
                colour = "black",
                size=4,
                vjust = 1.5
              )  +annotate(
                "text",
                x = -5,
                y = -21.5,
                label = "x1",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -4,
                y = -21.5,
                label = "x2",
                colour = "black",
                size=4,
                vjust = 1.5
              ) +annotate(
                "text",
                x = -3,
                y = -21.5,
                label = "x3",
                colour = "black",
                size=4,
                vjust = 1.5
              ) + annotate(
                "text",
                x =-2,
                y =  -21.5,
                label = "x4",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =-1,
                y =  -21.5,
                label = "x5",
                colour = "black",
                size=4,
                vjust = 1.5
              )+ annotate(
                "text",
                x =0,
                y =  -21.5,
                label = "x6",
                colour = "black",
                size=4,
                vjust = 1.5
              )+
              labs(x = "X", y = "F(X)") +
              theme_bw()+
              geom_line(data = xIntegracionComp2df) +
              geom_line(color = "blue", size = 1) 
            
            output$graficaIntegracionComp2 <- renderPlot(grafComp2)
            
          }
          if(input$botonPMComp2==5){
            output$AnaliticoIntegracionComp2<-renderText(paste( withMathJax("$$\\int_{-6}^{0}x^3+6x^2dx = \\frac{(x^3·(x+8))}{4} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
                                     $$\\int_{-6}^{0} f(x)dx = F(0)-F(-6) =\\frac{(0^3·(0+8)}{4} - \\frac{(-6)^3·(-6+8)}{4} = 108 $$ ")))
            
          }
      })
      
      
}

shinyApp(ui = ui, server = server)

