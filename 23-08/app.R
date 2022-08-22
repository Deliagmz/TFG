
library (shiny)
library(shinythemes)
library(ggplot2)
library(SciViews)
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
  
  return (result)
}


RPuntoMedioC <-function(F, a, b,N) {
  h <-(b-a)/(N+2)

  sum<-0
  for(j in 0: (N/2)){
    sum<- sum+ feval(F, a+2*j*h)
   
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

polinomioInt<- function(fila1,columna1){
  result<-""
 
  for(fil in 2:dim(fila1)[2]){
  
    if(fil==2){
      result<-paste(result,fila1[fil], sep = " ")
    }else if(fil!=1) {
      result<-paste(result,"+(",fila1[fil],sep = " ")
     for(numF in 3:fil){
        result<-paste(result,"(x -",columna1[numF-2],")", sep = " ")
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
      result<-paste(result,fila1[fil], sep = " ")
    }else if(fil!=1) {
      result<-paste(result,"+(",fila1[fil],sep = " ")
      for(numF in 3:fil){
        result<-paste(result,"*(x -",columna1[numF-2],")", sep = " ")
      }
      result<-paste(result,")",sep = " ")
    }
    
    
  }
  return (result)
}


intNewton_table <- function(x, y) {
  
  c <- length(x)
  f<-length(x)
  
  z <- pracma::zeros(c,  c+1)
  
  z[, 1] <- t(x)
  z[, 2] <- t(y)
 
  for (col in 3:(c+1)) {
    
    for (fil in 1:(f-1)) {
      
        z[fil, col] <- (z[fil + 1, col - 1] - z[fil,
                                              col - 1])/(z[fil + col - 2] - z[fil]) 
        
       
        if( z[fil, col]==Inf ||  z[fil, col]==-Inf|| is.nan(z[fil, col] )){
            z[fil, col]=0
           
        }
    }
   
    f<-f-1
    
  }
  print(z)
  return(z)
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
                   "Esta web tiene como proposito ayudar a estudiantes de forma interactiva a resolver
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
      tabPanel("Polinómica", uiOutput("IPolinomica"),
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
        uiOutput("INewton"),
        tabsetPanel(
          tabPanel("Teoría", uiOutput("TINewton"),
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
                                     &&&& f[x_0,x_1]\\\\
      1 & x_1 & f(x_1)  & f[x_1]&& f[x_0,x_1,x_2]\\\\
                                     &&&& f[x_1,x_2]&& f[x_0,x_1,x_2,x_3]\\\\
      2 & x_2& f(x_2)  & f[x_2]&& f[x_1,x_2,x_3]\\\\
                                     &&&& f[x_2,x_3]\\\\
      3 & x_3& f(x_3)  & f[x_3]
      \\end{array}$$"
                                   )
                                 ),
                                 p("Aunque en algunos casos también se calcula de la siguiente manera:"),
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
                                   "En la primera fila de la tabla se observan los coeficientes del polinomio de interpolación"
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
              h4("Interpolación de Newton en diferencias divididas"),
              br(),
            ),
            fluidRow(
              splitLayout(
                cellWidths = 500,
                DT::dataTableOutput("IntNew"),
                #tableOutput("IntNew"),
                plotOutput("graficaIntNew")
              ),
              br(),
              fluidRow(
                splitLayout(
                  cellWidths = 500,
                  fluidPage(
                    verticalLayout(
                      p("Introduce los puntos y pulsa siguiente paso cuando hayas terminado : "),
                      numericInput("x", "X:", 0, min = -50, max = 50),
                      numericInput("y", "Y:", 0, min = -50, max = 50),
                      splitLayout(
                        actionButton("botonIntNewP", "Introducir punto"),
                        shinyjs::useShinyjs(),
                        id = "botonInterpol",
                        actionButton("botonIntNew", "Sig. Paso", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")))),
                  htmlOutput("textoIntNew")
                ),
                
                
                
              ),
              
            )),
          ))
      ),
      tabPanel(
        "Inversa",
        uiOutput("IInversa"),
        tabsetPanel(
          tabPanel("Teoría", uiOutput("TIInversa"),
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
              h4("Interpolación inversa"),
              br(),
            ),
            fluidRow(
              splitLayout(
                cellWidths = 500,
                DT::dataTableOutput("IntInv"),
                plotOutput("graficaIntInv")
              ),
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
                        actionButton("botonIntInvP", "Introducir punto"),
                        actionButton("botonIntInv", "Sig. Paso", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")),
                      p("Introduce el valor para resolver la ecuación : "),
                      numericInput("y2", "Valor:", 0, min = -50, max = 50), 
                      actionButton("botonIntInvV", "Introducir valor") )),
                  
                  fluidPage(
                    verticalLayout(    
                      htmlOutput("textoIntInv"),
                      htmlOutput("textoIntInvValor")))
                ),
                
                
                
              ),
              
            )),
          ))
      ),
      tabPanel(
        "Osculatoria",
        uiOutput("IOculatoria"),
        tabsetPanel(
          tabPanel("Teoría", uiOutput("TIOsculatoria"),
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
          tabPanel("Ejemplos", uiOutput("EIOsculatoria")),
        )
      )
    ),
    
    
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
        uiOutput("ICerradaNewton"),
        tabsetPanel(
          tabPanel("Teoría", uiOutput("teoríaICN"),
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
            uiOutput("ECNewton"),
            mainPanel(navlistPanel(
              
              tabPanel("Ejemplo 1 ",mainPanel(
                
                fluidRow(
                  splitLayout(cellWidths = 325,h4(""),
                              h4("  $$\\int_{2}^{4}\\frac{x-2}{\\sqrt{(x-1)}}dx  $$ "), h4("")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasNCC", "Selecciona un método:",
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
                        uiOutput("textNCC"),
                        actionButton("botonNCC", "",icon = icon("sync"), style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                      )),
                    plotOutput("graficaNCC",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
              tabPanel("Ejemplo 2",mainPanel(
                
                fluidRow(
                  splitLayout(cellWidths = 325,h4(""),
                              h4(" $$\\int_{0}^{4}x\\sqrt{(x^2+9)}dx $$ "), h4("")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasNCS", "Selecciona un método:",
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
                        uiOutput("textNCS"),
                        uiOutput("AnaliticoNCS"),
                        actionButton("botonNCS", "",icon = icon("sync"), style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                      )),
                    plotOutput("graficaNCS",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
              
            ))
          )
        )),
      
      
      
      tabPanel(
        "Fórmula abierta de Newton-Côtes",
        uiOutput("IAbiertaNewton"), tabsetPanel(
          tabPanel("Teoría", uiOutput("teoríaIAN"),
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
            uiOutput("EANewton"),
            mainPanel(navlistPanel(
              tabPanel("Ejemplo 1 ",mainPanel(
                
                fluidRow(
                  splitLayout(cellWidths = 325,h4(""),
                              h4(" $$\\int_{3}^{5}\\frac{x^3}{(x^2+1)(x^2+1)(x^2+9)}dx $$ "), h4("")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasNAC", "Selecciona un método:",
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
                        uiOutput("textNAC"),
                        actionButton("botonNAC", "",icon = icon("sync"), style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                      )),
                    plotOutput("graficaNAC",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
              tabPanel("Ejemplo 2 ",mainPanel(
                
                fluidRow(
                  splitLayout(cellWidths = 325,h4(""),
                              h4("  $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx $$  "), h4("")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasNAS", "Selecciona un método:",
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
                        uiOutput("DatosNAS"),
                        uiOutput("TablaNAS"),
                        uiOutput("textNAS"),
                        uiOutput("resultNAS"),
                        uiOutput("AnaliticoNAS"),
                       
                          uiOutput("botonPMNAS"),
                          uiOutput("boton2PNAS"),
                          uiOutput("boton3PNAS")
                        
                      )),
                    plotOutput("graficaNAS",height = 600,width = 500),
                    
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
          tabPanel("Teoría", uiOutput("teoríaIC"),
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
                
                fluidRow(
                  splitLayout(cellWidths = 400,
                            h4(" $$\\int_{0}^{2}|x-2|^3(1-sen(\\pi x))dx $$ "), h4("Considerando los subintervalos: $$[0,\\frac{1}{2}],
                                                                                   [\\frac{1}{2},1],[1,\\frac{3}{2}], [\\frac{3}{2},2]$$")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasCC", "Selecciona un método de integración compuesta:",
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
                        uiOutput("DatosCC"),
                        uiOutput("TablaCC"),
                        uiOutput("textCC"),
                        uiOutput("resultCC"),
                          uiOutput("botonTrapecioCC"),
                          uiOutput("botonSimpsonCC"),
                          uiOutput("botonPMCC")
                        
                      )),
                    plotOutput("graficaCC",height = 600,width = 500),
                    
                  ),
                  
                  
                ),
                
              )),
              tabPanel("Ejemplo 2 ",mainPanel(
                
                fluidRow(
                  splitLayout(cellWidths = 325,h4(""),
                              h4("  $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx $$  "), h4("")),
                  br(),
                  splitLayout(cellWidths = 700,
                              fluidPage( 
                                verticalLayout(
                                  radioButtons(
                                    "ReglasCS", "Selecciona un método:",
                                    c("Regla del trapecio compuesta", "Regla de Simpson compuesta", "Regla del punto medio compuesta"),
                                    selected = "Regla del trapecio compuesta",
                                    inline = TRUE
                                  ),
                                  
                                )),
                              
                              
                              
                  ),br()),
                fluidRow(
                  splitLayout(
                    cellWidths = 600,
                    fluidPage( 
                      verticalLayout(
                        uiOutput("textCS"),
                        uiOutput("AnaliticoCS"),
                        actionButton("botonCS", "",icon = icon("sync"), style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                      )),
                    plotOutput("graficaCS",height = 600,width = 500),
                    
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
            round(dfNewtonExp2[i, 1], 2),
            "-\\frac{",
            round(dfNewtonExp2[i, 2], 2) ,
            "}{",
            round(dfNewtonExp2[i, 3], 2),
            "} = ",
            round(dfNewtonExp2[i, 4], 2) ,
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
            round(dfNewtonTrig2[i, 1], 4),
            "-\\frac{",
            round(dfNewtonTrig2[i, 2], 2) ,
            "}{",
            round(dfNewtonTrig2[i, 3], 2),
            "} = ",
            round(dfNewtonTrig2[i, 4], 2) ,
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
            round(dfNewtonPC2[i, 1], 4),
            "-\\frac{",
            round(dfNewtonPC2[i, 2], 2) ,
            "}{",
            round(dfNewtonPC2[i, 3], 2),
            "} = ",
            round(dfNewtonPC2[i, 4], 2) ,
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
            round(dfNewtonPC2[i, 1], 4),
            "-\\frac{",
            round(dfNewtonPC2[i, 2], 2) ,
            "}{",
            round(dfNewtonPC2[i, 3], 2),
            "} = ",
            round(dfNewtonPC2[i, 4], 2) ,
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
            round(dfNewtonOpt2[i, 1], 4),
            "-\\frac{",
            round(dfNewtonOpt2[i, 2], 2) ,
            "}{",
            round(dfNewtonOpt2[i, 3], 2),
            "} = ",
            round(dfNewtonOpt2[i, 4], 2),
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
  
  
  
  ##INTERPOLACION NEWTON
  
  ##INTRODUCIR DATOS
  observe( v$m <- data.frame(x = numeric(), "f(x)" = numeric()))
  cont<-0
  observeEvent(input$botonIntNewP,{
    
    req(input$x,input$y) 
    
    tmp <- data.frame(x = input$x,y  = input$y)
    colnames (tmp) <- c ( 'X', 'F(X)')
    
    size<-dim(v$m)
    if(size[2]<3){
      output$textoIntNew<-renderText(paste(""))
      
      v$m <- rbind(v$m,tmp)
      output$IntNew <- DT::renderDataTable({
        v$m
      }) 
      gInterN<-ggplot((v$m)) +
        geom_point(aes(x =v$m[,1],y=v$m[,2]))+theme_bw()+ labs(x = "X", y = "F(X)")
      output$graficaIntNew <-
        renderPlot(gInterN)
    }
  })
  
  ##REALIZAR TABLA
  
  
  counter <- reactiveValues(countervalue = 0)
  ##ELIMINAR DATOS QUE SOBRAN
  observeEvent(input$botonIntNew,{
    counter$countervalue <- counter$countervalue + 1   
    tabla<- as.data.frame(v$m)
    
    tam<-dim(tabla)
    
    if(tam[2]<tam[1]+1){
      tab<- intNewton_table(tabla[,1],tabla[,2])
      tam2<-dim(tab)  
      if(tab[,3:tam2[2]]== pracma::zeros(tam2[1],tam2[2]-2)){ 
        v$m <- data.frame(x = numeric(), "f(x)" = numeric())
        gInterInf<-ggplot( v$m )
        output$graficaIntNew <-
          renderPlot(gInterInf)
        output$textoIntNew<-renderText(paste("<h4 style=color:#FF0000;>ERROR <h4> <h4>Datos no validos, por favor introduce nuevos puntos. <h4>"))
        counter$countervalue <-0  
      }else{
        fx_n<- tab[, counter$countervalue+2]
        v$m <- cbind(v$m,
                     fx_n)
        colnames(v$m)[counter$countervalue+2]<-paste("fx_",counter$countervalue+2)
        
        output$IntNew <- DT::renderDataTable({
          v$m()%>% formatRound(columns=c('X', 'F(X)'), digits=3)
        })
        
      }
    }else{
      polinomio<-polinomioInt(v$m[1,],v$m[,1])
      output$textoIntNew<-renderText(paste("<h4>El polinomio de interpolación es: <br><br> P(x) =", polinomio," <h4>"))
      dataX <- as.data.frame(v$m[,1])
      if(dim(v$m)[1]<3){
        gInterInf<-ggplot(v$m )+geom_line(aes(x=v$m[,1],y=v$m[,2]),color="blue",size=1)+ theme_bw()+ labs(x = "X", y = "F(X)")
        output$graficaIntNew <-
          renderPlot(gInterInf)
      }else{
        
        gInterInf<-ggplot(v$m)+geom_point(aes(x=v$m[,1],y=v$m[,2]),color="blue",size=1)+ theme_bw()+ labs(x = "X", y = "F(X)")
        output$graficaIntNew <-
          renderPlot(curva)
      }
      
    }
    
    
    
  })
  
  
  ##INTERPOLACION INVERSA
  
  ##INTRODUCIR DATOS
  observe( inv$m  <- data.frame(x = numeric(), "f(x)" = numeric()))
  observeEvent(input$botonIntInvP,{
    
    req(input$x1,input$y1) 
    
    tmpInv <- data.frame(x = input$x1,y  = input$y1)
    colnames (tmpInv) <- c ( 'X', 'F(X)')
    
    sizeInv<-dim(inv$m )
    if(sizeInv[2]<3){
      output$textoIntInv<-renderText(paste(""))
      
      inv$m <- rbind(inv$m ,tmpInv)
      output$IntInv <- DT::renderDataTable({
        inv$m 
      }) 
      gInterInv<-ggplot((inv$m )) +
        geom_point(aes(x =inv$m [,1],y=inv$m [,2]))+theme_bw()+ labs(x = "X", y = "F(X)")
      output$graficaIntInv <-
        renderPlot(gInterInv)
    }
  })
  
  ##REALIZAR TABLA
  
  
  counterInv <- reactiveValues(countervalueInv = 0,countervalueInv2 = 0)
  ##ELIMINAR DATOS QUE SOBRAN
  observeEvent(input$botonIntInv,{
    counterInv$countervalueInv <- counterInv$countervalueInv + 1   
    tablaInv<- as.data.frame(inv$m )
    
    tamInv<-dim(tablaInv)
    
    if(tamInv[2]<tamInv[1]+1){
      tabInv<- intNewton_table(tablaInv[,1],tablaInv[,2])
      
      if(tabInv[,3:tamInv[2]]== pracma::zeros(tamInv[1],tamInv[2])){ 
        inv$m  <- data.frame(x = numeric(), "f(x)" = numeric())
        gInterInvInf<-ggplot(inv$m )
        output$graficaIntInv <-
          renderPlot(gInterInvInf)
        output$textoIntInv<-renderText(paste("<h4 style=color:#FF0000;>ERROR <h4> <h4>Datos no validos, por favor introduce nuevos puntos. <h4>"))
        counterInv$countervalueInv <-0  
      }else{
        fx<- tabInv[, counterInv$countervalueInv+2]
        inv$m  <- cbind(inv$m ,
                        fx)
        colnames(inv$m)[counterInv$countervalueInv+2]<-paste("fx_",counter$countervalue+2)
        
        output$IntInv <- DT::renderDataTable({
          inv$m 
        })
        
      }
    }else{
      polinomioInv<-polinomioInt(inv$m [1,],inv$m [,1])
      output$textoIntInv<-renderText(paste("<h4>El polinomio de interpolación es: <br><br> P(x) =", polinomioInv," <h4>"))
      if(dim(inv$m )[1]<3){
        gInterInvInf<-ggplot(inv$m)+geom_line(aes(x=inv$m [,1],y=inv$m [,2]),color="blue",size=1)+ theme_bw()+ labs(x = "X", y = "F(X)")
        output$graficaIntInv <-
          renderPlot(gInterInvInf)
      }
      
    }
    
    
    
  })
  #BOTÓN VALOR
  observeEvent(input$botonIntInvV,{
    counterInv$countervalueInv2 <- counterInv$countervalueInv2 + 1   
    req(input$y2)
    y2<-input$y2
    if(counterInv$countervalueInv2<2){
      output$textoIntInvValor<-renderText(paste("El valor introducido es ",y2))
      tablaInv2<- as.data.frame(inv$m )
      tamInv2<-dim(tablaInv2)
      polinomioInv2<-polinomioInt(inv$m [1,],inv$m [,1])
      if(tamInv2[2]>=tamInv2[1]+1){
        output$textoIntInvValor<-renderText(paste("Resolvemos ", polinomioInv2," =",y2))
        
      }
    }
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
      output$textoDerivacion1 <- renderText( paste( withMathJax("$$f'(1.3) = \\frac{f(1.3+",h,")-f(1.3-",h,")}{2·",h,"} =$$$$\\frac{ f(",1.3+h,")-f(",1.3-h,")}{",2*h,"}= \\frac{ ",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3+h),2),
                                                        "-",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3-h),2),"}{",2*h,"}=",derC,"$$")))
      
    }else if(input$metodosDerivacion1== "Descentrada izquierda" ){
      output$textoDerivacion1 <- renderText( paste( withMathJax("$$f'(1.3) = \\frac{f(1.3)-f(1.3-",h,")}{",h,"} =$$$$\\frac{ f(",1.3,")-f(",1.3-h,")}{",h,"}= \\frac{ ",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3),2),
                                                        "-",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3-h),2),"}{",h,"}=",derI,"$$")))
      
    } else {
      output$textoDerivacion1 <-renderText( paste( withMathJax("$$f'(1.3) = \\frac{f(1.3+",h,")-f(1.3)}{",h,"} =$$$$\\frac{ f(",1.3+h,")-f(",1.3,")}{",h,"}= \\frac{ ",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3+h),2),
                                                       "-",round(feval(function(x) x^(exp(1)-x^2)*exp((-x)^2)*(-2*x*ln(x)+(1/x)),1.3),2),"}{",h,"}=",derD,"$$")))
      
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
      output$textoDerivacion2 <- renderText( paste( withMathJax("$$f'(1.3) = \\frac{f(1.3+",h2,")-f(1.3-",h2,")}{2·",h2,"} =$$ $$\\frac{ f(",1.3+h2,")-f(",1.3-h2,")}{",2*h2,"}= \\frac{ ",round(feval(function(x)(x^2+3*x-2)^4,1.3+h2),2),"-",round(feval(function(x)(x^2+3*x-2)^4,1.3-h2),2),"}{",2*h2,"}=",derC2,"$$")))
      
    }else if(input$metodosDerivacion2== "Descentrada izquierda" ){
      output$textoDerivacion2 <- renderText( paste( withMathJax("$$f'(1.3) = \\frac{f(1.3)-f(1.3-",h2,")}{",h2,"} =$$$$\\frac{ f(",1.3,")-f(",1.3-h2,")}{",h2,"}= \\frac{ ",round(feval(function(x)(x^2+3*x-2)^4,1.3),2),"-",round(feval(function(x)(x^2+3*x-2)^4,1.3-h2),2),"}{",h2,"}=",derI2,"$$")))
      
    } else {
      output$textoDerivacion2 <-renderText( paste( withMathJax("$$f'(1.3) = \\frac{f(1.3+",h2,")-f(1.3)}{",h2,"} =$$$$\\frac{ f(",1.3+h2,")-f(",1.3,")}{",h2,"}= \\frac{ ",round(feval(function(x)(x^2+3*x-2)^4,1.3+h2),2),"-",round(feval(function(x)(x^2+3*x-2)^4,1.3),2),"}{",h2,"}=",derD2,"$$")))
      
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
  
  
  ##INTEGRACIÓN NEWTONCOTES cerrada C
  
  observeEvent(input$botonNCC,{
    
    trapecioNCC<- Rtrapecio(function(x) 
      (x-2)/(x-1)^(1-2) ,2, 4)
    simpsonNCC<- Rsimpson(function(x)
      (x-2)/(x-1)^(1-2) ,2, 4)
    simpson3octNCC<- Rsimpson3oct(function(x)
      (x-2)/(x-1)^(1-2) ,2, 4)
    
    if(input$ReglasNCC== "Regla del trapecio" ){
      
      output$textNCC <- renderText( paste( withMathJax("Sabiendo que \\(x_0 = a = 2, x_n = b = 4\\) y \\(n= 1\\)  $$h =\\frac{b-a}{n}= \\frac{4 -2}{1}=2$$
                                     $$\\int_{2}^{4}\\frac{x-2}{\\sqrt(x-1)}dx  =\\frac{2}{2}(f(2)+f(4))= ",trapecioNCC,"$$")))
      
    }else if(input$ReglasNCC== "Regla de Simpson" ){
      output$textNCC <- renderText( paste( withMathJax("Sabiendo que \\( x_i=x_0+ih, x_0 = a = 2, x_n = b = 4\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n}= \\frac{4 -2}{2}=1$$
                                      $$x_1=2+1·1=3$$
                                     $$\\int_{2}^{4}\\frac{x-2}{\\sqrt(x-1)}dx  =\\frac{1}{3}(f(2)+4f(3)+f(4))= ",simpsonNCC,"$$")))
      
    } else {
      output$textNCC <-renderText( paste( withMathJax("Sabiendo que \\( x_i=x_0+ih, x_0 = a = 2, x_n = b = 4\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n}= \\frac{4 -2}{3}=\\frac{2}{3}$$
                                                       $$x_1=2+1·\\frac{2}{3}=\\frac{8}{3}$$  $$x_2=2+2·\\frac{2}{3}=\\frac{10}{3}$$
                                                      $$\\int_{2}^{4}\\frac{x-2}{\\sqrt(x-1)}dx  =\\frac{3·2}{8·3}(f(2)+3f(\\frac{8}{3})+3f(\\frac{10}{3})+f(4))= ",simpson3octNCC,"$$")))
      
    }
    
    xNCC <- seq(1.5, 5, 0.05)
    xNCCdf <- as.data.frame(xNCC)
    
    grafTrapC <- ggplot(xNCCdf, aes(x = xNCC, y =  (xNCC-2)/(xNCC-1)^(1-2)))+
      geom_ribbon(data = xNCCdf, aes(x = xNCC, xmax =4, xmin = 2), fill = "#f2540c",alpha=0.8)+
      geom_ribbon(data = xNCCdf, aes(x = xNCC, ymax =Inf , ymin =  (xNCC-2)/(xNCC-1)^(1-2)), fill = "white",alpha=0.7)+
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
      )+labs(x = "X", y = "F(X)") +
      theme_bw()+
      geom_line(data = xNCCdf) +
      geom_line(color = "blue", size = 1) 
    
    grafSimpsonC <- ggplot(xNCCdf, aes(x = xNCC, y =   (xNCC-2)/(xNCC-1)^(1-2)))+
      geom_ribbon(data = xNCCdf, aes(x = xNCC, xmax =4 , xmin = 2), fill = "#f2540c",alpha=0.8)+
      geom_ribbon(data = xNCCdf, aes(x = xNCC, ymax =Inf , ymin =  (xNCC-2)/(xNCC-1)^(1-2)), fill = "white",alpha=0.7)+
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
      )+labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xNCCdf) +
      geom_line(color = "blue", size = 1)
    
    grafSimpson3C <- ggplot(xNCCdf, aes(x = xNCC, y =   (xNCC-2)/(xNCC-1)^(1-2)))+
      geom_ribbon(data = xNCCdf, aes(x = xNCC, xmax =4 , xmin = 2), fill = "#f2540c",alpha=0.8)+
      geom_ribbon(data = xNCCdf, aes(x = xNCC, ymax =Inf , ymin =  (xNCC-2)/(xNCC-1)^(1-2)), fill = "white",alpha=0.7)+
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
      )+labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xNCCdf) +
      geom_line(color = "blue", size = 1)
    
    if(input$ReglasNCC== "Regla del trapecio" ){
      output$graficaNCC <- renderPlot(grafTrapC)
    }else if(input$ReglasNCC== "Regla de Simpson" ){
      output$graficaNCC <- renderPlot(grafSimpsonC)
    } else {
      output$graficaNCC <- renderPlot(grafSimpson3C)
    }
    
  }
  )
  
  ##INTEGRACIÓN NEWTONCOTES cerrada S
  observeEvent(input$botonNCS,{
    
    trapecio<- Rtrapecio(function(x) 
      x*(x^2+9)^(1/2) ,0, 4 )
    simpson<- Rsimpson(function(x)
      x*(x^2+9)^(1/2),0, 4 )
    simpson3oct<- Rsimpson3oct(function(x)
      x*(x^2+9)^(1/2) ,0, 4 )
    
    if(input$ReglasNCS== "Regla del trapecio" ){
      output$textNCS <- renderText( paste( withMathJax("Sabiendo que \\(x_0 = a = 0, x_n = b = 4\\) y \\(n= 1\\)  $$h =\\frac{b-a}{n}= \\frac{4-0}{1}$$
                                     $$\\int_{0}^{4} x\\sqrt{(x^2+9)} dx =\\frac{4}{2}(f(0)+f(4))= ",trapecio,"$$")))
      
    }else if(input$ReglasNCS== "Regla de Simpson" ){
      output$textNCS <- renderText( paste( withMathJax("Sabiendo que \\( x_i=x_0+ih, x_0 = a = 0, x_n = b = 4\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n}= \\frac{4-0}{2}=2$$
                                    $$x_1 = 0+1·2=2$$
                                     $$\\int_{0}^{4}x\\sqrt{(x^2+9)} dx =\\frac{2}{3}(f(0)+4f(2)+f(4))= ",simpson,"$$")))
      
    } else {
      output$textNCS <-renderText( paste( withMathJax("Sabiendo que \\( x_i=x_0+ih, x_0 = a = 0, x_n = b = 4\\) y \\(n= 3\\)  $$h =\\frac{b-a}{n}= \\frac{4-0}{3}$$
                                                        $$x_1 = 0+1·\\frac{4}{3}=\\frac{4}{3}$$ $$x_2=0+2·\\frac{4}{3}=\\frac{8}{3}$$
                                                      $$\\int_{0}^{4} x\\sqrt{(x^2+9)} dx =\\frac{3·4}{8·3}(f(0)+3f(\\frac{4}{3})+3f(\\frac{8}{3})+f(4))= ",simpson3oct,"$$")))
      
    }
    
    output$AnaliticoNCS<-renderText(paste( withMathJax("Resolución de forma analítica:  $$\\int_{0}^{4} x\\sqrt{(x^2+9)} dx = \\frac{(x^2+9)^{\\frac{3}{2}}}{3} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
                                     $$\\int_{0}^{4} f(x)dx = F(4)-F(0) =\\frac{(4^2+9)^{\\frac{3}{2}}}{3} - \\frac{(0^2+9)^{\\frac{3}{2}}}{3} = 32.67 $$ ")))
    
    
    xNCS <- seq(-1, 5, 0.05)
    xNCSdf <- as.data.frame(xNCS)
    
    grafTrapS <- ggplot(xNCSdf, aes(x = xNCS, y =  (xNCS*(xNCS^2+9)^(1/2))))+
      geom_ribbon(data = xNCSdf, aes(x = xNCS, xmax =4 , xmin = 0), fill = "#f2540c",alpha=0.8)+
      geom_ribbon(data = xNCSdf, aes(x = xNCS, ymax =Inf , ymin = (xNCS*(xNCS^2+9)^(1/2))), fill = "white",alpha=0.7)+ 
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
      )+labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xNCSdf) +
      geom_line(color = "blue", size = 1) 
    

    
    grafSimpsonS <- ggplot(xNCSdf, aes(x = xNCS, y =  (xNCS*(xNCS^2+9)^(1/2))))+
      geom_ribbon(data = xNCSdf, aes(x = xNCS, xmax =4 , xmin = 0), fill = "#f2540c",alpha=0.8)+
      geom_ribbon(data = xNCSdf, aes(x = xNCS, ymax =Inf , ymin = (xNCS*(xNCS^2+9)^(1/2))), fill = "white",alpha=0.7) +
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
      )+labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xNCSdf) +
      geom_line(color = "blue", size = 1)
    
    grafSimpson3S <- ggplot(xNCSdf, aes(x = xNCS, y =  (xNCS*(xNCS^2+9)^(1/2))))+
      geom_ribbon(data = xNCSdf, aes(x = xNCS, xmax =4 , xmin = 0), fill = "#f2540c",alpha=0.8)+
      geom_ribbon(data = xNCSdf, aes(x = xNCS, ymax =Inf , ymin = (xNCS*(xNCS^2+9)^(1/2))), fill = "white",alpha=0.7)+
      geom_ribbon(aes( ymax =0 , ymin = -Inf), fill = "white",alpha=0.7)+
      
      geom_line(aes(y=0,colour = "y"),color="black", size = 1,linetype = 1)+
      geom_line(aes(x = 0,colour = "x0"),color="black", size = 1,linetype = 5) +
      geom_line(aes(x = 1.34,colour = "x1"),color="black", size = 1,linetype = 5) +
      geom_line(aes(x = 2.7,colour = "x2"),color="black", size = 1,linetype = 5) +
      geom_line(aes(x = 4,colour = "x3"),color="black", size = 1,linetype = 5)+
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
      )+labs(x = "X", y = "F(X)") +theme_bw()+
      geom_line(data = xNCSdf) +
      geom_line(color = "blue", size = 1)  
    
    if(input$ReglasNCS== "Regla del trapecio" ){
      output$graficaNCS <- renderPlot(grafTrapS)
    }else if(input$ReglasNCS== "Regla de Simpson" ){
      output$graficaNCS <- renderPlot(grafSimpsonS)
    } else {
      output$graficaNCS <- renderPlot(grafSimpson3S)
    }
    
  }
  )
  

    ##INTEGRACIÓN NEWTONCOTES Abierta C
    
    observeEvent(input$botonNAC,{
      
      PuntoMedioNAC<- RPuntoMedio(function(x) 
        (x^3)/((x^2+1)*(x^2+4)*(x^2+9)) ,3, 5)
      abierta2pNAC<- fabierta2P(function(x)
        (x^3)/((x^2+1)*(x^2+4)*(x^2+9)) ,3, 5)
      abierta3pNAC<- fabierta3P(function(x)
        (x^3)/((x^2+1)*(x^2+4)*(x^2+9)) ,3, 5)
      
      if(input$ReglasNAC== "Regla del punto medio" ){
        output$textNAC <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a=3, x_{n+1} = b=5\\) y \\(n= 0\\)  $$h =\\frac{b-a}{n+2}= \\frac{5 -3}{2}=1$$
                                    $$x_0= a+h = 4$$
                                     $$\\int_{3}^{5}\\frac{x^3}{(x^2+1)(x^2+4)(x^2+9)}dx  =2·1(f(4))= ", PuntoMedioNAC,"$$")))
        
      }else if(input$ReglasNAC== "Fórmula de 2 puntos" ){
        output$textNAC <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a=3, x_{n+1} = b=5\\) y \\(n= 1\\)  $$h =\\frac{b-a}{n+2}= \\frac{5 -3}{1+2}=\\frac{2}{3}$$
                                                           $$x_0= a+h = 3+\\frac{2}{3}=\\frac{11}{3}$$
                                                           $$x_1= b-h = 5-\\frac{2}{3}=\\frac{13}{3}$$
                                                         $$\\int_{3}^{5}\\frac{x^3}{(x^2+1)(x^2+4)(x^2+9)}dx  = \\frac{3\\frac{2}{3}}{2}(f(\\frac{11}{3})+f(\\frac{13}{3}))= ",abierta2pNAC,"$$")))
        
      } else {
        output$textNAC <-renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a=3, x_{n+1} = b=5\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n+2}= \\frac{5 -3}{2+2}=\\frac{1}{2}$$
                                                         $$x_0= a+h = 3+\\frac{1}{2}=\\frac{7}{2}$$
                                                         $$x_1= x_0 +ih = \\frac{7}{2}+1·\\frac{1}{2}=4$$
                                                        $$x_2= b-h = 5-\\frac{1}{2}=\\frac{9}{2}$$
                                     $$\\int_{3}^{5}\\frac{x^3}{(x^2+1)(x^2+4)(x^2+9)}dx  =\\frac{4\\frac{1}{2}}{3}(2f(\\frac{7}{2})-f(4)+2f(\\frac{9}{2}))= ",abierta3pNAC,"$$")))
        
      }
      
      xNAC <- seq(0, 6, 0.05)
      xNACdf <- as.data.frame(xNAC)
      
      grafPuntoMC <- ggplot(xNACdf, aes(x = xNAC, y = (xNAC^3)/((xNAC^2+1)*(xNAC^2+4)*(xNAC^2+9))))+
        geom_ribbon(data = xNACdf, aes(x = xNAC , xmax =5 , xmin = 3), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xNACdf, aes(x = xNAC, ymax =Inf , ymin = ((xNAC^3)/((xNAC^2+1)*(xNAC^2+4)*(xNAC^2+9)))), fill = "white",alpha=0.7)+
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
        )+
        labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xNACdf) +
        geom_line(color = "blue", size = 1) 
        
      
      graf2PC<-ggplot(xNACdf, aes(x = xNAC, y = (xNAC^3)/((xNAC^2+1)*(xNAC^2+4)*(xNAC^2+9))))+
        geom_ribbon(data = xNACdf, aes(x = xNAC , xmax =5 , xmin = 3), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xNACdf, aes(x = xNAC, ymax =Inf , ymin = ((xNAC^3)/((xNAC^2+1)*(xNAC^2+4)*(xNAC^2+9)))), fill = "white",alpha=0.7)+
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
        )+labs(x = "X", y = "F(X)") +
        theme_bw()+
        geom_line(data = xNACdf) +
        geom_line(color = "blue", size = 1) 
      
      graf3PC<-ggplot(xNACdf, aes(x = xNAC, y = (xNAC^3)/((xNAC^2+1)*(xNAC^2+4)*(xNAC^2+9))))+
        geom_ribbon(data = xNACdf, aes(x = xNAC , xmax =5 , xmin = 3), fill = "#f2540c",alpha=0.8)+
        geom_ribbon(data = xNACdf, aes(x = xNAC, ymax =Inf , ymin = ((xNAC^3)/((xNAC^2+1)*(xNAC^2+4)*(xNAC^2+9)))), fill = "white",alpha=0.7)+
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
        geom_line(data = xNACdf) +
        geom_line(color = "blue", size = 1) 
   
      if(input$ReglasNAC== "Regla del punto medio" ){
        output$graficaNAC <- renderPlot(grafPuntoMC)
      }else if(input$ReglasNAC== "Fórmula de 2 puntos" ){
        output$graficaNAC <- renderPlot(graf2PC)
      } else {
        output$graficaNAC <- renderPlot(graf3PC)
      }
      
    }
    )
    
    ##INTEGRACIÓN NEWTONCOTES abierta S
    observeEvent(input$ReglasNAS,{
      
      output$DatosNAS<-NULL
      output$TablaNAS<-NULL
      output$textNAS<-NULL
      output$resultNAS<-NULL
      output$AnaliticoNAS<-NULL
      
      output$botonPMNAS<-NULL
      output$boton2PNAS<-NULL
      output$boton3PNAS<-NULL
      output$graficaNAS<-NULL
      
      if(input$ReglasNAS== "Regla del punto medio"){
        output$botonPMNAS <- renderUI({
          actionButton("botonPMNAS", label = ">>" ,style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
        
      }else if(input$ReglasNAS== "Fórmula de 2 puntos" ){
        output$boton2PNAS <- renderUI({
          actionButton("boton2PNAS", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
        
      } else {
        output$boton3PNAS <- renderUI({
          actionButton("boton3PNAS", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
      }
      
    })
    observeEvent(input$botonPMNAS,{
      PuntoMedioNAS<- RPuntoMedio(function(x) 
        cos(x)*sin(x) ,0, pi/2)
      xNAS <- seq(-0.25, 1.75, 0.05)
      xNASdf <- as.data.frame(xNAS)
      
      if(input$botonPMNAS==1){
        output$DatosNAS <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a=0, x_{n+1} = b=\\frac{\\pi}{2}\\) y \\(n= 0\\)  $$h =\\frac{b-a}{n+2}= \\frac{\\frac{\\pi}{2}-0}{2}=\\frac{\\pi}{4}$$")))
        grafPuntoMS<-ggplot(xNASdf, aes(x = xNAS, y = cos(xNAS)*sin(xNAS)))+
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
          geom_line(data = xNASdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaNAS <- renderPlot(grafPuntoMS)
         }
      if(input$botonPMNAS==2){
        output$TablaNAS <- renderText( paste( withMathJax("$$x_0= a+h =\\frac{\\pi}{4}$$ $$f(x_0)= cos(\\frac{\\pi}{4})sin(\\frac{\\pi}{4})=\\frac{1}{2}$$")))
        
        grafPuntoMS <- ggplot(xNASdf, aes(x = xNAS, y = cos(xNAS)*sin(xNAS)))+
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
          )+ labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xNASdf) +
          geom_line(color = "blue", size = 1) 
        
        
        output$graficaNAS <- renderPlot(grafPuntoMS)
      }
      if(input$botonPMNAS==3){
        output$textNAS <- renderText( paste( withMathJax("  $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx  =2·\\frac{\\pi}{4}(f(\\frac{\\pi}{4}))= $$ ")))
      }
      if(input$botonPMNAS==4){
        output$resultNAS <- renderText( paste( withMathJax("$$\\frac{\\pi}{2}·\\frac{1}{2}=", PuntoMedioNAS,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
        grafPuntoMS <- ggplot(xNASdf, aes(x = xNAS, y = cos(xNAS)*sin(xNAS)))+
          geom_ribbon(data = xNASdf, aes(x = xNAS, xmax =pi/2 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xNASdf, aes(x = xNAS, ymax =Inf , ymin = cos(xNAS)*sin(xNAS)), fill = "white",alpha=0.7)+
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
          )+ labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xNASdf) +
          geom_line(color = "blue", size = 1) 
        
        
        output$graficaNAS <- renderPlot(grafPuntoMS)
        
      }
      if(input$botonPMNAS==5){
      output$AnaliticoNAS<-renderText(paste( withMathJax(" $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx = - \\frac{cos(x)^2}{2} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
    $$\\int_{0}^{\\frac{\\pi}{2}} f(x)dx = F(\\frac{\\pi}{2})-F(0) = - \\frac{cos(\\frac{\\pi}{2})^2}{2} - (-\\frac{cos(0)^2}{2})= 0.5 $$ ")))
      }
    })
    
    observeEvent(input$boton2PNAS,{
      abierta2pNAS<- fabierta2P(function(x)
        cos(x)*sin(x) ,0, pi/2)
      xNAS <- seq(-0.25, 1.75, 0.05)
      xNASdf <- as.data.frame(xNAS)
      
      if(input$boton2PNAS==1){
        output$DatosNAS <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a, x_{n+1} = b , x_0= a+h ,x_n= b-h \\) y \\(n= 1\\)  $$h =\\frac{b-a}{n+2}= \\frac{\\frac{\\pi}{2}-0}{1+2}=\\frac{\\pi}{6}$$ ")))
        graf2PS<-ggplot(xNASdf, aes(x = xNAS, y = cos(xNAS)*sin(xNAS)))+
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
          geom_line(data = xNASdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaNAS <- renderPlot(graf2PS)
        }
      if(input$boton2PNAS==2){
        output$TablaNAS <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
       i & -1 & 0 & 1 & 2   \\\\
       \\hline
       x_i & 0 & \\frac{\\pi}{6} &\\frac{\\pi}{3} & \\frac{\\pi}{2}\\\\
     
      f(x_i) & 0  & \\frac{\\sqrt{3}}{4} & \\frac{\\sqrt{3}}{4} & 0 \\\\
                                   
      \\end{array}$$")))
        
        graf2PS<-ggplot(xNASdf, aes(x = xNAS, y =cos(xNAS)*sin(xNAS)))+
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
          )+labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xNASdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaNAS <- renderPlot(graf2PS)
        
      }
      if(input$boton2PNAS==3){
        output$textNAS <- renderText( paste( withMathJax("$$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx = \\frac{3\\frac{\\pi}{6}}{2}(f(\\frac{\\pi}{6})+f(\\frac{\\pi}{3}))= $$ ")))
      }
      if(input$boton2PNAS==4){
        output$resultNAS <- renderText( paste( withMathJax("$$\\frac{\\pi}{4}(\\frac{\\sqrt{3}}{4}+\\frac{\\sqrt{3}}{4})=", abierta2pNAS,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
        graf2PS<-ggplot(xNASdf, aes(x = xNAS, y =cos(xNAS)*sin(xNAS)))+
          geom_ribbon(data = xNASdf, aes(x = xNAS, xmax =pi/2 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xNASdf, aes(x = xNAS, ymax =Inf , ymin = cos(xNAS)*sin(xNAS)), fill = "white",alpha=0.7)+
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
          )+labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xNASdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaNAS <- renderPlot(graf2PS)
         
      }
      if(input$boton2PNAS==5){
        output$AnaliticoNAS<-renderText(paste( withMathJax(" $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx = - \\frac{cos(x)^2}{2} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
    $$\\int_{0}^{\\frac{\\pi}{2}} f(x)dx = F(\\frac{\\pi}{2})-F(0) = - \\frac{cos(\\frac{\\pi}{2})^2}{2} - (-\\frac{cos(0)^2}{2})= 0.5 $$ ")))
      }
    })
    
    
    observeEvent(input$boton3PNAS,{
      abierta3pNAS<- fabierta3P(function(x)
        cos(x)*sin(x) ,0, pi/2)
      xNAS <- seq(-0.25, 1.75, 0.05)
      xNASdf <- as.data.frame(xNAS)
      
      if(input$boton3PNAS==1){
        output$DatosNAS <- renderText( paste( withMathJax("Sabiendo que \\(x_{-1} = a, x_{n+1} = b,x_0= a+h,x_i= x_0+ih, x_n= b-h\\) y \\(n= 2\\)  $$h =\\frac{b-a}{n+2}= \\frac{\\frac{\\pi}{2}-0}{2+2}=\\frac{\\pi}{8}$$")))
        graf3PS<-ggplot(xNASdf, aes(x = xNAS, y = cos(xNAS)*sin(xNAS)))+
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
          geom_line(data = xNASdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaNAS <- renderPlot(graf3PS)
      }
      if(input$boton3PNAS==2){
        output$TablaNAS <- renderText( paste( withMathJax(" $$\\begin{array}{|c|c|c|c|c|}
                                                                        i & -1 & 0 & 1 & 2 & 3  \\\\
                                                                      \\hline
                                                                      x_i & 0 & \\frac{\\pi}{8} &\\frac{\\pi}{4} &\\frac{3\\pi}{8} & \\frac{\\pi}{2}\\\\
                                                                      
                                                                      f(x_i) & 0  & 0.354  & \\frac{\\pi}{2} & 0.354 &0 \\\\
                                                                      
                                                                      \\end{array}$$")))
        
        
        graf3PS<-ggplot(xNASdf, aes(x = xNAS, y = cos(xNAS)*sin(xNAS)))+
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
          )+labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xNASdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaNAS <- renderPlot(graf3PS)
        
      }
      if(input$boton3PNAS==3){
        output$textNAS <- renderText( paste( withMathJax("$$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx =\\frac{4\\frac{\\pi}{8}}{3}(2f(\\frac{\\pi}{8})-f(\\frac{\\pi}{4})+2f(\\frac{3\\pi}{8}))=$$")))
      }
      if(input$boton3PNAS==4){
        output$resultNAS <- renderText( paste( withMathJax("$$\\frac{\\pi}{6}(2·0.354-\\frac{1}{2}+2·0.354)=", abierta3pNAS,"$$ $$ \\textbf{Resolver de forma analítica:}$$")))
        
        graf3PS<-ggplot(xNASdf, aes(x = xNAS, y = cos(xNAS)*sin(xNAS)))+
          geom_ribbon(data = xNASdf, aes(x = xNAS, xmax =pi/2 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xNASdf, aes(x = xNAS, ymax =Inf , ymin = cos(xNAS)*sin(xNAS)), fill = "white",alpha=0.7)+
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
          )+labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xNASdf) +
          geom_line(color = "blue", size = 1) 
        output$graficaNAS <- renderPlot(graf3PS)
       
      }
      if(input$boton3PNAS==5){
        output$AnaliticoNAS<-renderText(paste( withMathJax(" $$\\int_{0}^{\\frac{\\pi}{2}} cos(x)sin(x)dx = - \\frac{cos(x)^2}{2} + C $$ Aplicando la regla de Barrow:  $$\\int_{a}^{b} f(x)dx = F(b)-F(a)$$
    $$\\int_{0}^{\\frac{\\pi}{2}} f(x)dx = F(\\frac{\\pi}{2})-F(0) = - \\frac{cos(\\frac{\\pi}{2})^2}{2} - (-\\frac{cos(0)^2}{2})= 0.5 $$ ")))
      }
    })
    
    
    
    
    
    ##INTEGRACIÓN COMPUESTA
    observeEvent(input$ReglasCC,{
    
      output$DatosCC<-NULL
      output$TablaCC<-NULL
      output$textCC<-NULL
      output$resultCC<-NULL
      
      output$botonTrapecioCC<-NULL
      output$botonSimpsonCC<-NULL
      output$botonPMCC<-NULL
      output$graficaCC<-NULL
      
      if(input$ReglasCC== "Regla del trapecio" ){
        output$botonTrapecioCC <- renderUI({
          actionButton("botonTrapecioCC", label = ">>" ,style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
        
      }else if(input$ReglasCC== "Regla de Simpson" ){
        output$botonSimpsonCC <- renderUI({
          actionButton("botonSimpsonCC", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
        
      } else {
        output$botonPMCC <- renderUI({
          actionButton("botonPMCC", label = ">>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
        })
      }
      
    })
    
    observeEvent(input$botonTrapecioCC,{
      trapecioComp<- RtrapecioC(function(x) 
        abs(x-2)^3*(1-sin(pi*x)),0, 2,4)
      
      xNCComp <- seq(-0.2, 2.5, 0.05)
      xNCCompdf <- as.data.frame(xNCComp)
      
      if(input$botonTrapecioCC==1){
        output$DatosCC <- renderText( paste( withMathJax("Sabiendo que \\(x_0=a=0, x_N=b=2\\) y \\(N= 4 \\)   $$h =\\frac{b-a}{N}= \\frac{2 -0}{4}=\\frac{1}{2}$$")))
        grafComp <- ggplot(xNCCompdf, aes(x = xNCComp, y =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))))+
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
          geom_line(data = xNCCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaCC <- renderPlot(grafComp)
      }
      if(input$botonTrapecioCC==2){
        output$TablaCC <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
       i & 0 & 1 &2 & 3 &4 \\\\
       \\hline
       x_i & 0 & 0.5 &1 & 1.5 & 2\\\\
     
      f(x_i) & 8  & 0  & 1 & \\frac{1}{4} &0 \\\\
                                   
      \\end{array}$$")))
        
        grafComp <- ggplot(xNCCompdf, aes(x = xNCComp, y =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))))+
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
            label = "a",
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
            label = "b",
            colour = "black",
            size=4,
            vjust = 1.5
          )+
          labs(x = "X", y = "F(X)") +
          theme_bw()+
          geom_line(data = xNCCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaCC <- renderPlot(grafComp)
        
      }
      if(input$botonTrapecioCC==3){
        output$textCC <- renderText( paste( withMathJax("$$\\int_{0}^{2}|x-2|^3(1-sen(\\pi x))dx =\\frac{\\frac{1}{2}}{2}(f(0)+2(f(0.5)+f(1)+f(1.5))+f(2))=$$ ")))
        
      }
      if(input$botonTrapecioCC==4){
        output$resultCC <- renderText( paste( withMathJax("$$\\frac{1}{4}(8+2(0+1+\\frac{1}{4})+0)=", trapecioComp,"$$")))
        grafComp <- ggplot(xNCCompdf, aes(x = xNCComp, y =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))))+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =0.5 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =1 , xmin = 0.5), fill = "green",alpha=0.8)+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =1.5 , xmin = 1), fill = "yellow",alpha=0.8)+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =2 , xmin = 1.5), fill = "purple",alpha=0.8)+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp, ymax =Inf , ymin =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))), fill = "white",alpha=0.7)+
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
          geom_line(data = xNCCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaCC <- renderPlot(grafComp)
        
      }
      })
      
      observeEvent(input$botonSimpsonCC,{
      simpsonComp<- RsimpsonC(function(x)
        abs(x-2)^3*(1-sin(pi*x)),0, 2,4)
      
      xNCComp <- seq(-0.2, 2.5, 0.05)
      xNCCompdf <- as.data.frame(xNCComp)
      if(input$botonSimpsonCC==1){
        output$DatosCC <- renderText( paste( withMathJax("Sabiendo que \\(x_0=a=0, x_N=b=2\\) y \\(N= 4\\)  $$h =\\frac{b-a}{N}= \\frac{2 -0}{4}=\\frac{1}{2}$$")))
        grafComp <- ggplot(xNCCompdf, aes(x = xNCComp, y =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))))+
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
          geom_line(data = xNCCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaCC <- renderPlot(grafComp)
      }
      if(input$botonSimpsonCC==2){
        output$TablaCC <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
                                                                 i & 0 & 1 &2 & 3 &4 \\\\
                                                                 \\hline
                                                                  x_i & 0 & \\frac{1}{2} &1 & \\frac{3}{2} & 2\\\\
                                                               
                                                                f(x_i) & 8  & 0  & 1 & \\frac{1}{4} &0 \\\\
                                                                                             
                                                                \\end{array}$$")))
        grafComp <- ggplot(xNCCompdf, aes(x = xNCComp, y =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))))+
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
          geom_line(data = xNCCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaCC <- renderPlot(grafComp)
        
      }
      if(input$botonSimpsonCC==3){
        output$textCC <- renderText( paste( withMathJax("$$\\int_{0}^{2}|x-2|^3(1-sen(\\pi x))dx =\\frac{\\frac{1}{2}}{3}(f(0)+2(f(1))+4(f(\\frac{1}{2})+f(\\frac{3}{2}))+f(2))= $$")))
      }
      
      if(input$botonSimpsonCC==4){
        output$resultCC <- renderText( paste( withMathJax("$$\\frac{1}{6}(8+2(1)+4(0+\\frac{1}{4})+0)= ",simpsonComp,"$$")))
        grafComp <- ggplot(xNCCompdf, aes(x = xNCComp, y =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))))+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =0.5 , xmin = 0), fill = "#f2540c",alpha=0.8)+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =1 , xmin = 0.5), fill = "green",alpha=0.8)+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =1.5 , xmin = 1), fill = "yellow",alpha=0.8)+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =2 , xmin = 1.5), fill = "purple",alpha=0.8)+
          geom_ribbon(data = xNCCompdf, aes(x = xNCComp, ymax =Inf , ymin =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))), fill = "white",alpha=0.7)+
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
          geom_line(data = xNCCompdf) +
          geom_line(color = "blue", size = 1) 
        
        output$graficaCC <- renderPlot(grafComp)
      }
      })
      
      observeEvent(input$botonPMCC,{
        puntoMedioComp<- RPuntoMedioC(function(x)
          abs(x-2)^3*(1-sin(pi*x)),0, 2,4)
        
        xNCComp <- seq(-0.2, 2.5, 0.05)
        xNCCompdf <- as.data.frame(xNCComp)
        if(input$botonPMCC==1){
          output$DatosCC <- renderText( paste( withMathJax("Sabiendo que \\(x_0=a=0, x_N=b=2\\) y \\(N= 4\\)  $$h =\\frac{b-a}{N+2}= \\frac{2 -0}{4+2}=\\frac{1}{3}$$")))
          grafComp <- ggplot(xNCCompdf, aes(x = xNCComp, y =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))))+
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
            geom_line(data = xNCCompdf) +
            geom_line(color = "blue", size = 1) 
          
          output$graficaCC <- renderPlot(grafComp)
          }
        if(input$botonPMCC==2){
          output$TablaCC <- renderText( paste( withMathJax("$$\\begin{array}{|c|c|c|c|c|}
                                                                 i & 0 & 1 &2 & 3 &4 \\\\
                                                                 \\hline
                                                                 x_i & 0 & \\frac{1}{2} &1 & \\frac{3}{2} & 2\\\\
                                                               
                                                                f(x_i) & 8  & 0  & 1 & \\frac{1}{4} &0 \\\\
                                                                                             
                                                                \\end{array}$$")))
          grafComp <- ggplot(xNCCompdf, aes(x = xNCComp, y =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))))+
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
            geom_line(data = xNCCompdf) +
            geom_line(color = "blue", size = 1) 
          
          output$graficaCC <- renderPlot(grafComp)
          
        }
        if(input$botonPMCC==3){
          output$textCC <- renderText( paste( withMathJax("$$\\int_{0}^{2}|x-2|^3(1-sen(\\pi x))dx =2·\\frac{1}{3}(f(0)+f(1)+f(2))= $$")))
        }
        
        if(input$botonPMCC==4){
          output$resultCC <- renderText( paste( withMathJax("$$\\frac{2}{3}(8+1+0)= ",puntoMedioComp,"$$")))
          grafComp <- ggplot(xNCCompdf, aes(x = xNCComp, y =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))))+
            geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =0.5 , xmin = 0), fill = "#f2540c",alpha=0.8)+
            geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =1 , xmin = 0.5), fill = "green",alpha=0.8)+
            geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =1.5 , xmin = 1), fill = "yellow",alpha=0.8)+
            geom_ribbon(data = xNCCompdf, aes(x = xNCComp , xmax =2 , xmin = 1.5), fill = "purple",alpha=0.8)+
            geom_ribbon(data = xNCCompdf, aes(x = xNCComp, ymax =Inf , ymin =  abs(xNCComp-2)^3*(1-sin(pi*xNCComp))), fill = "white",alpha=0.7)+
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
            geom_line(data = xNCCompdf) +
            geom_line(color = "blue", size = 1) 
          
          output$graficaCC <- renderPlot(grafComp)
          
        }
      })
}

shinyApp(ui = ui, server = server)

