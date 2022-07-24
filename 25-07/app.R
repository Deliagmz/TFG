
library (shiny)
library(shinythemes)
library(ggplot2)
Sys.setlocale(category = "LC_ALL", locale = "Spanish")
##FUNCIONES BIPART
feval <- function(f, ...) do.call(f, list(...))

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
  res <- data.frame(x0 = NULL,
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
                   x_n = x ,
                   `f(x)` = Fx,
                   `f'(x)`=DFx,
                   `x_n+1` = x1
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
    tabPanel(
      "MÉTODOS NUMÉRICOS",
      uiOutput("opcion0"),
      mainPanel(fluidRow(
        column(
          8,
          br(),
          h4(
            "Esta web tiene como proposito ayudar a estudiantes de forma interactiva a resolver
         y entender diferentes problemas relacionados con la asignatura MÉTODOS NUMÉRICOS. "
          ),
          br(),
          h4(
            "Puedes consultar información y ejercicios resueltos sobre el temario que se encuentra en el menu superior."
          ),
          h4(
            "Los ejercicios se resolveran paso a paso al ritmo que tu impongas, tanto los cálculos como la gráfica que se genera a partir de estos."
          ),
          h4(
            "Encontraras tanto ejemplos ya creados como otros donde seras tú el que deberá introducir los datos del ejercicio que deseas resolver."
          ),
          
          h3("¿Comenzamos?"),
          br(),
        ),
        column(1,
               offset = 1,
               br(),
               br(),
               br(),
               div(
                 style = "display: block;",
                 img(
                   src = "descarga.jpg",
                   height = 200,
                   width = 250,
                   align = "center"
                 )
               )),
        br(),
        br()  ,
      ),
      
      fluidRow(column(
        12,
        br(),
        h5("Trabajo fin de grado realizado por Delia Gómez Lobato."),
        h5("Grado en Ingeniería Informática."),
        br(),
        img(
          src = "uma2.png",
          height = 100,
          width = 125,
          align = "center"
        ),
        img(
          src = "rstudio.png",
          height = 35,
          width = 100,
          align = "center"
        )
        
      )))
    ),
    ##APARTADO ECUACIONES LINEALES
    tabPanel(
      "Ecuaciones no lineales",
      tabsetPanel(
        tabPanel("Teoría", uiOutput("teoria"), mainPanel(fluidRow
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
                                                               p("1) Busco un intervalo \\([a , b ]\\) tal que \\(f(a)f(b)<0\\)"),
                                                               br(),
                                                               p("2) Repetir hasta criterio de parada:"),
                                                               column(
                                                                 12,
                                                                 offset = 0.5,
                                                                 p("Tomo  \\(x=\\frac{a+b}{2}\\)."),
                                                                 p(
                                                                   "Evalúo \\(f(x):\\begin{cases}
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
                                                             p("Partimos de un punto inicial \\(x_0\\) e iteramos mediante: "),
                                                             column(
                                                               12,
                                                               offset = 2,
                                                               p("\\(x_k+1= x_k -\\frac{f(x_k)}{f'(x_k)} \\)", style =
                                                                   "font-size:large"),
                                                             ),
                                                           ),
                                                         ))),
        tabPanel(
          "Ejemplos Bipartición",
          uiOutput("Biparticion"),
          mainPanel(navlistPanel(
            tabPanel("Ecuación exponencial",  mainPanel(
              fluidRow(
                br(),
                h4("Ecuación exponencial con intervalo [0,1]: $$f(x)=e^{(3x)}-4$$ "),
                br(),
              ),
              fluidRow(
                splitLayout(
                  cellWidths = 500,
                  fluidPage(
                    verticalLayout(
                      tableOutput("exp"), htmlOutput("textexp"))),
                  
                  plotOutput("graficaExp",height = 600,width = 500)
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botonexp", ">>>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
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
                  fluidPage(
                    verticalLayout(
                      tableOutput("trig"), htmlOutput("texttrig"))),
                  
                  plotOutput("graficaTrig",height = 600,width = 500)
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botontrig", ">>>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                ),
                
              ),
              
            )),
            
            tabPanel("Punto de corte entre funciones",  mainPanel(
              fluidRow(
                br(),
                
                h4(
                  "Punto de corte entre las funciones $$f(x)=x^2-4x+5$$ $$g(x)=x+1$$ "
                ),
                h4(
                  " $$f(x)=g(x)$$ $$F(x)=f(x)-g(x)=0$$ "
                ),
                h4(
                  "Ecuación con intervalo [2,7]: $$F(x)=x^2-5x+4$$ "
                ),
                br(),
              ),
              fluidRow(
                splitLayout(
                  cellWidths = 500,fluidPage( verticalLayout(
                    tableOutput("ptoscorte"), htmlOutput("textPC"))),
                  plotOutput("graficaPC",height = 600,width = 500)
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botonPC", ">>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                ),
                
              ),
              
            )),
            
            
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
                  cellWidths = 500, fluidPage(
                    verticalLayout(
                      tableOutput("opt"), htmlOutput("textOpt"))),
                  plotOutput("graficaOpt",height = 600,width = 500)
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botonopt", ">>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                ),
                
              ),
              
            )),
            
          ))
        ),
        tabPanel(
          "Ejemplos Newton",
          uiOutput("Newton"),
          mainPanel(navlistPanel(
            tabPanel("Ecuación exponencial",  mainPanel(
              fluidRow(
                br(),
                splitLayout(cellWidths = 500,
                h4("Ecuación exponencial  $$f(x)=e^{(x)}+x^2-4$$ $$f'(x)=e^{(x)}+2x$$ con \\(x_0=3.5\\)"), 
                fluidPage( 
                  verticalLayout(br(),br(),p("Recordamos:"), 
                                 p("\\(x_k+1= x_k -\\frac{f(x_k)}{f'(x_k)} \\)", style = "font-size:large"))),
                
              ),br()),
              fluidRow(
                splitLayout(
                 cellWidths = 500,
                 fluidPage(
                 verticalLayout(
                  tableOutput("expN"), htmlOutput("textExpN"))),
                 
                  plotOutput("graficaExpN",height = 600,width = 500),
                
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botonexpN", ">>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                ),
                
              ),
              
            )),
            
            tabPanel("Ecuación trigonométrica", mainPanel(
              fluidRow(
                br(),
                splitLayout(cellWidths = 500,
                            h4("Ecuación trigonométrica $$f(x)=x^2 -cos(x)$$$$f'(x)=2x+sen(x)$$ con \\(x_0=2\\) "), 
                            fluidPage( 
                              verticalLayout(br(),br(),p("Recordamos:"), 
                                             p("\\(x_k+1= x_k -\\frac{f(x_k)}{f'(x_k)} \\)", style = "font-size:large"))),
                            
                ),br()),
              fluidRow(
                splitLayout(
                  cellWidths = 500,
                  fluidPage(
                    verticalLayout(
                      tableOutput("trigN"), htmlOutput("textTrigN"))),
                  
                  plotOutput("graficaTrigN",height = 600,width = 500),
                  
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botonTrigN", ">>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                ),
                
              ),
              
            )),
            
            tabPanel("Punto de corte entre funciones",  mainPanel(
              fluidRow(
                br(),
                splitLayout(cellWidths = 500,
                h4(
                  "Punto de corte entre las funciones $$f(x)=x^4+3x^2-10$$ $$g(x)=sen(x)-x^3$$  $$f(x)=g(x)$$ $$F(x)=f(x)-g(x)=0$$"
                ),
               fluidPage( 
                verticalLayout(br(),br(), h4(
                  "Con \\(x_0=-0.5\\) $$F(x)=x^4-x^3+3x^2-sen(x)-10$$ $$F'(x)=4x^3-3x^2+6x-cos(x)$$ "
                ),p("Recordamos:"), 
                               p("\\(x_k+1= x_k -\\frac{f(x_k)}{f'(x_k)} \\)", style = "font-size:large"))),
              
           ), br()),
              fluidRow(
                splitLayout(
                  cellWidths = 500,fluidPage( verticalLayout(
                    tableOutput("ptoscorteN"), htmlOutput("textPCN"))),
                  plotOutput("graficaPCN",height = 600,width = 500)
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botonPCN", ">>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                ),
                
              ),
              
            )),
            
            tabPanel("Optimización",  mainPanel(
              fluidRow(
                br(),
                splitLayout(cellWidths = 500,
                            h4("Ejercicio optimización $$f(x)=x^3+4x^2-10$$ $$f'(x)=3x^2+8x$$ $$f''(x)=6x+8$$ con \\(x_0=0.75\\) "), 
                            fluidPage( 
                              verticalLayout(br(),br(),p("En este caso:"), 
                                             p("\\(x_k+1= x_k -\\frac{f'(x_k)}{f''(x_k)} \\)", style = "font-size:large"))),
                            
                ),br()),
              fluidRow(
                splitLayout(
                  cellWidths = 500,
                  fluidPage(
                    verticalLayout(
                      tableOutput("optN"), htmlOutput("textOptN"))),
                  
                  plotOutput("graficaOptN",height = 600,width = 500),
                  
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botonOptN", ">>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                ),
                
              ),
              
            )),
            
          )),
        )
      )
    ),
    
    ##APARTADO INTERPOLACIÃN
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
          tabPanel("Ejemplos", uiOutput("EIInversa")),
        ),
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
                                 strong("Interpolacion de Hermite", style = "color: #f2540c"),
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
               tabPanel("Teoría", uiOutput("teoría"),
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
                 uiOutput("EDerivacion"),
                 mainPanel(navlistPanel(
                   tabPanel("Fórmula centrada", "cetrada"),
                   tabPanel("Fórmula descentrada derecha", "izq"),
                   tabPanel("Fórmula descentrada izquierda", "drch"),
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
              tabPanel("Regla del trapecio", "cetrada"),
              tabPanel("Regla de Simpson", "izq"),
              tabPanel("Regla de Simpson tres octavos", "drch"),
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
                         "\\(x_i=x_0+ih\\) , \\( i = 0,...,n \\) con \\(x_0=a+h\\) , \\( x_n =b-h(x_{-1}=a,x_{n+1}=b)\\) , y \\(h = \\frac{b-a}{n+2} \\)"),
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
              tabPanel(" Regla del punto medio", "cetrada"),
              tabPanel("Fórmula abierta de dos puntos", "izq"),
              tabPanel("Fórmula abierta de tres puntos", "drch"),
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
                         "Siendo \\(h=\\frac{b-a}{N+2}, N =2m\\), teniendo en cuenta que  \\((\\frac{h^2(b-a)}{6}f''(\\xi) )\\) representa el error y que (\\( x_{-1}<\\xi< x_3\\)):"
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
              tabPanel("Regla del trapecio compuesta", "cetrada"),
              tabPanel("Regla de Simpson compuesta", "izq"),
              tabPanel("Regla del punto medio compuesta", "drch"),
            ))
          )
          
        )
        
      ),
      
    )
  )
  
)
  v <- reactiveValues(m=data.frame(x = numeric(), "f(x)" = numeric()))
 v$m<-NULL
 v$m <- data.frame(x = numeric(), "f(x)" = numeric())

server <- function(input, output) { 

  ##BIPARTICION
  #EXPONENCIAL
 
  bipartExp <- function(i, n,ab) {
    if (i <= n) {
      output$exp <- renderTable(df[1:i, ],digits=3)
      output$graficaExp <-
        renderPlot(
          gBipartEXP +  
            geom_line(aes(x = ab[i, 1],colour = "a"),color="#f2540c", size = 1) +
            geom_line(aes(x = ab[i, 2],colour = "b"),color="#f2540c", size = 1)
          
          +annotate(
              "text",
              x = ab[i, 1]-0.05,
              y = 15,
              label = "a",
              colour = "black",
              vjust = 1.5
            )
          + annotate(
              "text",
              x =ab[i, 2]+0.05,
              y =  15,
              label = "b",
              colour = "black",
              vjust = 1.5
            )
        )
    }
    if(i==n){
      output$textexp<-renderText(paste("<h4> La raíz de la ecuación se encuentra entre <br>[0.4609375, 0.46484375] con cota de error 0.004 </h4>"))
      output$graficaExp <-
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
              x = 0.46+0.05,
              y = 0,
              label ="Raíz",
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
  }
  x<-seq(0,1,0.05)
  df <- bipart_table(function(x)
    (exp(3 * x) - 4), 0, 1, 8, 1)
  colnames (df) <- c ('n', 'X_n', 'F (X_n)','[a, b]', 'Cota de error')
  ab <- bipart(function(x)
    (exp(3 * x) - 4), 0, 1, 8)

 xDT <- as.data.frame(x)
  gBipartEXP <- ggplot(xDT, aes(x = x, y =  exp(3 * x) - 4)) +
    labs(x = "X", y = "F(X)") +
    theme(axis.title.x = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +
    theme(axis.title.y = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +theme_bw()+
    geom_line(data = xDT) +
    geom_line(color = "blue", size = 1)
  
 
  output$graficaExp <- renderPlot(gBipartEXP)
  
  
  observeEvent(input$botonexp, bipartExp(input$botonexp, 8, ab))
  
  #TRIGONOMETRICA
  
  bipartTrig <- function(i, n,abT) {
    if (i <= n) {
      output$trig <- renderTable(dfT[1:i, ],digits=3)
      output$graficaTrig <-
        renderPlot(
          gBipartTRIG +  
            geom_line(aes(x = abT[i, 1],colour = "a"),color="#f2540c", size = 1) +
            geom_line(aes(x = abT[i, 2],colour = "b"),color="#f2540c", size = 1)
          
          +annotate(
            "text",
            x = abT[i, 1]-0.1,
            y = 1.3,
            label = "a",
            colour = "black",
            vjust = 1.5
          )
          + annotate(
            "text",
            x =abT[i, 2]+0.1,
            y =  1.3,
            label = "b",
            colour = "black",
            vjust = 1.5
          )
        )
    }
    if(i==n){
      output$texttrig<-renderText(paste("<h4> La raíz de la ecuación se encuentra entre <br>[-0.8515625, -0.84375] con cota de error 0.008 </h4>"))
      output$graficaTrig <-
        renderPlot(
          gBipartTRIG 
          + annotate(
            "point",
            x =-0.849 ,
            y = 0,
            colour = "red",
            size = 2
          ) +
            annotate(
              "text",
              x = -0.849+0.1,
              y = 0,
              label ="Raíz",
              colour = "black",
              vjust = 1.5
            )
        )
    }
  }
  
  xT<-seq(-1,1,0.05)
  dfT <- bipart_table(function(xT)
    sin(xT) + cos(xT ^ 2), -1, 1, 8, 1)
  colnames (dfT) <- c ('n', 'X_n', 'F (X_n)','[a, b]', 'Cota de error')
   abT <- bipart(function(xT)
    sin(xT) + cos(xT ^ 2), -1, 1, 8)
  
 
  
  xDTrig <- as.data.frame(xT)
  gBipartTRIG <- ggplot(xDTrig, aes(x = xT, y =  sin(xT) + cos(xT ^ 2))) +
    labs(x = "X", y = "F(X)") +
    theme(axis.title.x = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +
    theme(axis.title.y = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +theme_bw()+
    geom_line(data = xDTrig) +
    geom_line(color = "blue", size = 1)
  
  output$graficaTrig <- renderPlot(gBipartTRIG)
  
  observeEvent(input$botontrig, bipartTrig(input$botontrig, 8, abT))
  
  
  #PUNTO DE CORTE 2 FUNCIONES
  
  bipartPC <- function(i, n,abPC) {
    if (i <= n) {
      output$ptoscorte <- renderTable(dfPC[1:i, ],digits=3)
      output$graficaPC <-
        renderPlot(
          gBipartPC +
            geom_line(aes(x = abPC[i, 1],colour = "a"),color="#f2540c", size = 1) +
            geom_line(aes(x = abPC[i, 2],colour = "b"),color="#f2540c", size = 1)
          
          +annotate(
            "text",
            x = abPC[i, 1]-0.15,
            y = 16,
            label = "a",
            colour = "black",
            vjust = 1.5
          )
          + annotate(
            "text",
            x =abPC[i, 2]+0.15,
            y =  16,
            label = "b",
            colour = "black",
            vjust = 1.5
          )
        )
    }
    if(i==n){
      output$textPC<-renderText(paste("<h4> El punto de corte entre f(x) y g(x) se encuentra en el <br>intervalo [3.9921875, 4.01171875] con cota de error 0.02 <br> <br>En este
                                       caso, en F(x)=0 -> x=4 <br> Al evaluar ese valor en g(x)=x+1 -> g(x)=5<br> <br>  El punto de corte es P(4,5) </h4>"))
      output$graficaPC <-
        renderPlot(
          gBipartPC +  
            annotate(
              "point",
              x =4 ,
              y = 5,
              colour = "blue",
              size = 2
            ) +
            annotate(
              "text",
              x = 4-0.8,
              y = 6.5,
              label ="Punto de corte",
              colour = "black",
              vjust = 1.5
            )
        )
    }
  }
  
  xPC<-seq(1.5,7,0.05)
  dfPC <- bipart_table(function(xPC)
    (xPC^2-5*xPC+4), 2, 7, 8, 1)
  colnames (dfPC) <- c ('n', 'X_n', 'F (X_n)','[a, b]', 'Cota de error')
   abPC <- bipart(function(xPC)
    (xPC^2-5*xPC+4), 2, 7, 8)
  
  
  
  xDPC <- as.data.frame(xPC)
  gBipartPC <- ggplot(xDPC, aes(x = xPC, y =   xPC^2-5*xPC+4)) + annotate(
    "text",
    x =7,
    y =  14.5,
    label = "F(x)",
    colour = "black",
    vjust = 1.5,size=5
  ) +geom_line( aes(x = xPC, y =   (xPC^2-4*xPC+5)),color="black",size=1)+ annotate(
    "text",
    x =7,
    y =  23,
    label = "f(x)",
    colour = "black",
    vjust = 1.5,size=5
  )+geom_line( aes(x = xPC, y =   (xPC+1)),color=" red",size=1)+ annotate(
      "text",
      x =7,
      y =  7,
      label = "g(x)",
      colour = "black",
      vjust = 1.5,size=5
    )+
    labs(x = "X", y = "F(X)") +
    theme(axis.title.x = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +
    theme(axis.title.y = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +theme_bw()+
    geom_line(data = xDPC) +
    geom_line(color = "blue", size = 1)
  
  output$graficaPC <- renderPlot(gBipartPC)
  
  observeEvent(input$botonPC, bipartPC(input$botonPC, 8, abPC))
  
  
  #OPTIMIZACION
  
  bipartOpt <- function(i, n,abO) {
    if (i <= n) {
      output$opt <- renderTable(dfO[1:i, ],digits=3)
      output$graficaOpt <-
        renderPlot(
          gBipartOPT +  
            geom_line(aes(x = abO[i, 1],colour = "a"),color="#f2540c", size = 1) +
            geom_line(aes(x = abO[i, 2],colour = "b"),color="#f2540c", size = 1)
          
          +annotate(
            "text",
            x = abO[i, 1]-0.1,
            y = 10,
            label = "a",
            colour = "black",
            vjust = 1.5
          )
          + annotate(
            "text",
            x =abO[i, 2]+0.1,
            y =  10,
            label = "b",
            colour = "black",
            vjust = 1.5
          )
        )
    }
    if(i==n){
      output$textOpt<-renderText(paste("<h4> El mínimo o máximo de la funcion f(x) se encuentra en el <br>intervalo [0.5703125, 0.578125] con cota de error 0.008 <br> <br>En este
                                       caso, en f'(x)=0 ->  x=0.57736 <br> Al evaluar ese valor en f''(x)=6x -> f''(x)=3.46<br> Se trata de un punto mínimo al ser un valor positivo.<br>
                                     <br>  El mínimo de la función f(x) es P(0.57736,-1.3849) </h4>"))
      output$graficaOpt <-
        renderPlot(
          gBipartOPT +  
            annotate(
              "point",
              x =0.57736 ,
              y = -1.3849,
              colour = "red",
              size = 2
            ) +
            annotate(
              "text",
              x = 0.57736+0.1,
              y = -1.3849,
              label ="Punto mínimo",
              colour = "black",
              vjust = 1.5
            )
        )
    }
  }

  xO<-seq(0,2,0.05)
  dfO <- bipart_table(function(xO)
    (3 * xO ^ 2 - 1), 0, 2, 8, 1)
  colnames (dfO) <- c ('n', 'X_n', 'F (X_n)','[a, b]', 'Cota de error')
  abO <- bipart(function(xO)
    (3 * xO ^ 2 - 1), 0, 2, 8)
  
  xDOpt <- as.data.frame(xO)
  gBipartOPT <- ggplot(xDOpt, aes(x = xO, y =   (3 * xO ^ 2 - 1)))+ annotate(
    "text",
    x =2,
    y =  4.5,
    label = "f(x)",
    colour = "black",
    vjust = 1.5,size=5
  ) +geom_line( aes(x = xO, y =   (xO ^ 3 -xO - 1)),color="black",size=1)+ annotate(
    "text",
    x =2,
    y =  10,
    label = "f '(x)",
    colour = "black",
    vjust = 1.5,size=5
  ) +
    labs(x = "X", y = "F(X)'") +
    theme(axis.title.x = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +
    theme(axis.title.y = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +theme_bw()+
    geom_line(data = xDOpt) +
    geom_line(color = "blue", size = 1)
  
  output$graficaOpt <- renderPlot(gBipartOPT)
  observeEvent(input$botonopt, bipartOpt(input$botonopt, 8, abO))
  
  
  
  ##NEWTON
  #EXPONENCIAL
  newtonExp <- function(i, n) {
    if (i <= n) {
      output$expN <- renderTable(dfNII[1:i,],digits=3)
      output$textExpN<-renderText(paste("<h4> x",i,"= ", round(dfNII[i,1],2), "-(", round(dfNII[i,2],2) ,"/", round(dfNII[i,3],2),") = ", round(dfNII[i,4],2) ,"<br>","</h4>"))
      output$graficaExpN <-
        renderPlot(
          g +   geom_point(aes(x = round(dfN[i + 1, 2],2)),color="black", size = 1) +
            
            geom_abline(
            slope = slope( round(dfN[i,2],2)),
            intercept = intercept(round(dfN[i,2],2)),
            color = "#f2540c",
            size = 1
          )
          + annotate(
            "point",
            x =round(dfN[i,2],2) ,
            y = round(dfN[i,4],2),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfN[i,2],2)-0.1,
              y =round(dfN[i,4],2),
              label =round(dfN[i,2],2),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfN[i + 1, 2],2) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfN[i + 1, 2],2)+0.15,
              y = 0,
              label =round( dfN[i + 1, 2],2),
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    if(i==n){
      output$textExpN<-renderText(paste("<h4> La raíz de la ecuación es 1.058 </h4>"))
      output$graficaExpN <-
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
              x = 1.058006+0.1,
              y = 0,
              label ="Raíz",
              colour = "black",
              vjust = 1.5
            )
        )
    }
    
    
  }
  
  dfNII <- newtonII(function(x)
    (exp(x) + x ^ 2 - 4),  function(x)
      (exp( x) + 2 * x)
    , 3.5)
  colnames (dfNII) <- c ( 'X_n', 'F(X_n)','F`(X_n)', 'X_(n+1)')
  dfN <- newton(function(x)
    (exp( x) + x ^ 2 - 4),  function(x)
      ( exp( x) + 2 * x)
    , 3.5)
  xN <- seq(0, 3.5, 0.05)
  
  
  #GRAFICA
  
  xNT <- as.data.frame(xN)
  g <- ggplot(xNT, aes(x = xN, y =  exp( xN) + xN ^ 2 - 4)) +
    labs(x = "X", y = "F(X)") +
    theme(axis.title.x = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +
    theme(axis.title.y = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +theme_bw()+
    geom_line(data = xNT) +
    geom_line(color = "blue", size = 1)
  
  
  #Y=mx+b  m=> SLOPE b=>INTERCEPT
  slope <- function(x) {
    xincrement <- 0.2
    yincrement <-
      (exp( (x + 0.1)) + (x + 0.1) ^ 2 - 4) - (exp( (x - 0.1)) + (x - 0.1) ^ 2 - 4)
    slope <- yincrement / xincrement
    return(slope)
  }
  
  intercept <- function(x) {
    intercept <- (exp( x) + x ^ 2 - 4) - slope(x) * x
    yvalues <- slope(x) * xN + intercept
    xNT <- as.data.frame(cbind(xNT, yvalues))
    return(intercept)
  }
  
  output$graficaExpN <- renderPlot(g)
  
  observeEvent(input$botonexpN, newtonExp(input$botonexpN, 6))
  
  
  #TRIGONOMÉTRICA
  newtonTrig <- function(i, n) {
    if (i <= n) {
      output$trigN <- renderTable(dfTrig2[1:i,],digits=4)
      output$textTrigN<-renderText(paste("<h4> x",i,"= ", round(dfTrig2[i,1],4), "-(", round(dfTrig2[i,2],2) ,"/", round(dfTrig2[i,3],2),") = ", round(dfTrig2[i,4],2) ,"<br>","</h4>"))
      output$graficaTrigN <-
        renderPlot(
          gTrig + geom_point(aes(x = round(dfTrig1[i + 1, 2],2)),color="black", size = 1) +
            geom_abline(
              slope = slopeTrig( dfTrig1[i,2]),
              intercept = interceptTrig(dfTrig1[i,2]),
              color = "#f2540c",
              size = 1
            )
          + annotate(
            "point",
            x =round(dfTrig1[i,2],4) ,
            y = round(dfTrig1[i,4],4),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfTrig1[i,2],4)-0.1,
              y =round(dfTrig1[i,4],4),
              label =round(dfTrig1[i,2],4),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfTrig1[i + 1, 2],4) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfTrig1[i + 1, 2],4)+0.15,
              y = 0,
              label =round( dfTrig1[i + 1, 2],4),
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    if(i==n){
      output$textTrigN<-renderText(paste("<h4> La raíz de la ecuación es 0.8241 </h4>"))
      output$graficaTrigN <-
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
              x = 0.824132+0.1,
              y = 0,
              label ="Raíz",
              colour = "black",
              vjust = 1.5
            )
        )
    }
    
  }
  
  dfTrig2 <- newtonII(function(x)
    (x^2-cos(x)),  function(x)
      (2*x+sin(x))
    , 2)
  colnames (dfTrig2) <- c ( 'X_n', 'F(X_n)','F`(X_n)', 'X_(n+1)')
  dfTrig1 <- newton(function(x)
    (x^2-cos(x)),  function(x)
      (2*x+sin(x))
    , 2)
  xNTrig <- seq(0, 2, 0.05)
  
  
  #GRAFICA
  
  xNTrigdf <- as.data.frame(xNTrig)
  gTrig <- ggplot(xNTrigdf, aes(x = xNTrig, y = (xNTrig^2-cos(xNTrig)))) +
    labs(x = "X", y = "F(X)") +
    theme(axis.title.x = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +
    theme(axis.title.y = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +theme_bw()+
    geom_line(data = xNTrigdf) +
    geom_line(color = "blue", size = 1)
  
  
  #Y=mx+b  m=> SLOPE b=>INTERCEPT
  slopeTrig <- function(x) {
    xincrementT <- 0.2
    yincrementT <-
      ((x+0.1)^2-cos(x+0.1))-((x-0.1)^2-cos(x-0.1))
    slopeT <- yincrementT / xincrementT
    return(slopeT)
  }
  
  interceptTrig <- function(x) {
    interceptT <- (x^2-cos(x)) - slopeTrig(x) * x
    yvaluesT <- slopeTrig(x) * xNTrig + interceptT
    xNTrigdf <- as.data.frame(cbind(xNTrigdf, yvaluesT))
    return(interceptT)
  }
  
  output$graficaTrigN <- renderPlot(gTrig)
  
  observeEvent(input$botonTrigN, newtonTrig(input$botonTrigN, 5))
  
  #PUNTO DE CORTE
  newtonPC <- function(i, n) {
    if (i < 4) {
      output$ptoscorteN <- renderTable(dfPC2[1:i,],digits=4)
      output$textPCN<-renderText(paste("<h4> x",i,"= ", round(dfPC2[i,1],4), "-(", round(dfPC2[i,2],2) ,"/", round(dfPC2[i,3],2),") = ", round(dfPC2[i,4],2) ,"<br>","</h4>"))
      output$graficaPCN <-
        renderPlot(
          gPC + geom_point(aes(x = round(dfPC1[i + 1, 2],2)),color="black", size = 1) +
            geom_abline(
              slope = slopePC( dfPC1[i,2]),
              intercept = interceptPC(dfPC1[i,2]),
              color = "#f2540c",
              size = 1
            )
          + annotate(
            "point",
            x =round(dfPC1[i,2],4) ,
            y = round(dfPC1[i,4],4),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfPC1[i,2],4)+0.2,
              y =round(dfPC1[i,4],4)+0.3,
              label =round(dfPC1[i,2],4),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfPC1[i + 1, 2],4) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfPC1[i + 1, 2],4)-0.2,
              y = -0.12,
              label =round( dfPC1[i + 1, 2],4),
              colour = "black",
              vjust = 1.5
            )
        )
      
    }else if(i>3 && i<n){
      output$ptoscorteN <- renderTable(dfPC2[1:i,],digits=4)
      output$textPCN<-renderText(paste("<h4> x",i,"= ", round(dfPC2[i,1],4), "-(", round(dfPC2[i,2],2) ,"/", round(dfPC2[i,3],2),") = ", round(dfPC2[i,4],2) ,"<br>","</h4>"))
      output$graficaPCN <-
        renderPlot(
          gPC + geom_point(aes(x = round(dfPC1[i + 1, 2],2)),color="black", size = 1) +
            geom_abline(
              slope = slopePC( dfPC1[i,2]),
              intercept = interceptPC(dfPC1[i,2]),
              color = "#f2540c",
              size = 1
            )
          + annotate(
            "point",
            x =round(dfPC1[i,2],4) ,
            y = round(dfPC1[i,4],4),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfPC1[i,2],4)-0.5,
              y =round(dfPC1[i,4],4)-0.7,
              label =round(dfPC1[i,2],4),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfPC1[i + 1, 2],4) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfPC1[i + 1, 2],4)+0.4,
              y = 3,
              label =round( dfPC1[i + 1, 2],4),
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    if(i==n){
      output$textPCN<-renderText(paste(" <h4>Hemos calculado que x = -1.5348<br> Evaluamos ese valor en la función g(x)-> g(x)=2.616 <br> El punto de corte entre f(x) y g(x) es P(-1.5348,2.616)</h4>"))
      output$graficaPCN <-
        renderPlot(
          gPC  +geom_line( aes(x = xNPC, y =   (xNPC^4+3*xNPC^2-10)),color="black",size=1)+ annotate(
            "text",
            x =-3.5,
            y =  150,
            label = "f(x)",
            colour = "black",
            vjust = 1.5,size=5
          )+geom_line( aes(x = xNPC, y =   (sin(xNPC)-xNPC^3)),color=" red",size=1)+ annotate(
            "text",
            x =-3.5,
            y =  30,
            label = "g(x)",
            colour = "black",
            vjust = 1.5,size=5
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
              x = -1.5348+0.5,
              y = 10,
              label ="Punto de corte",
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    
    
  }
  
  dfPC2 <- newtonII(function(x)
    (x^4+x^3+3*x^2-sin(x)-10),  function(x)
      (4*x^3+3*x^2+6*x-cos(x))
    , -0.5)
  colnames (dfPC2) <- c ( 'X_n', 'F(X_n)','F`(X_n)', 'X_(n+1)')
  dfPC1 <- newton(function(x)
    (x^4+x^3+3*x^2-sin(x)-10),  function(x)
      (4*x^3+3*x^2+6*x-cos(x))
    , -0.5)
  xNPC <- seq(-3.5, 1, 0.05)
  
  
  #GRAFICA
  
  xNPCdf <- as.data.frame(xNPC)
  gPC <- ggplot(xNPCdf, aes(x = xNPC, y =xNPC^4+xNPC^3+3*xNPC^2-sin(xNPC)-10 )) + annotate(
    "text",
    x =-3.5,
    y =  110,
    label = "F(x)",
    colour = "black",
    vjust = 1.5,size=5
  )+
    labs(x = "X", y = "F(X)") +
    theme(axis.title.x = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +
    theme(axis.title.y = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +theme_bw()+
    geom_line(data = xNPCdf) +
    geom_line(color = "blue", size = 1)
  
  
  #Y=mx+b  m=> SLOPE b=>INTERCEPT
  slopePC <- function(x) {
    xincrementPC <- 0.2
    yincrementPC <-
      ((x+0.1)^4+(x+0.1)^3+3*(x+0.1)^2-sin((x+0.1))-10)-((x-0.1)^4+(x-0.1)^3+3*(x-0.1)^2-sin((x-0.1))-10)
    slopepto <- yincrementPC / xincrementPC
    return(slopepto)
  }
  
  interceptPC <- function(x) {
    interceptpto <- (x^4+x^3+3*x^2-sin(x)-10) - slopePC(x) * x
    yvaluesPC <- slopePC(x) * xNPC + interceptpto
    xNPCdf <- as.data.frame(cbind(xNPCdf, yvaluesPC))
    return(interceptpto)
  }
  
  output$graficaPCN <- renderPlot(gPC)
  
  observeEvent(input$botonPCN, newtonPC(input$botonPCN, 7))
  
  
  
  #OPTIMIZACIÓN
  newtonOpt <- function(i, n) {
    if (i <= n) {
      output$optN <- renderTable(dfOpt2[1:i,],digits=4)
      output$textOptN<-renderText(paste("<h4> x",i,"= ", round(dfOpt2[i,1],4), "-(", round(dfOpt2[i,2],2) ,"/", round(dfOpt2[i,3],2),") = ", round(dfOpt2[i,4],2) ,"<br>","</h4>"))
      output$graficaOptN <-
        renderPlot(
          gOpt + geom_point(aes(x = round(dfOpt1[i + 1, 2],2)),color="black", size = 1) +
            geom_abline(
              slope = slopeOpt( dfOpt1[i,2]),
              intercept = interceptOpt(dfOpt1[i,2]),
              color = "#f2540c",
              size = 1
            )
          + annotate(
            "point",
            x =round(dfOpt1[i,2],4) ,
            y = round(dfOpt1[i,4],4),
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfOpt1[i,2],4)-0.1,
              y =round(dfOpt1[i,4],4),
              label =round(dfOpt1[i,2],4),
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = round(dfOpt1[i + 1, 2],4) ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = round(dfOpt1[i + 1, 2],4)+0.15,
              y = 0,
              label =round( dfOpt1[i + 1, 2],4),
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    if(i==n){
        output$textOptN<-renderText(paste("<h4> x5=0 -> f''(0)=6*0+8=8 <br> <br>El resultado positivo implica un punto mínimo <br><br> El mínimo de la funcion f(x) es P(0,-10)  </h4>"))
        output$graficaOptN <-
          renderPlot(
            gOpt  +geom_abline(
              slope = slopeOpt( dfOpt1[i,2]),
              intercept = interceptOpt(dfOpt1[i,2]),
              color = "#f2540c",
              size = 1
            )
            + annotate(
              "point",
              x =round(dfOpt1[i,2],4) ,
              y = round(dfOpt1[i,4],4),
              colour = "black",
              size = 2
            ) +
              annotate(
                "text",
                x = round(dfOpt1[i,2],4)-0.1,
                y =round(dfOpt1[i,4],4),
                label =round(dfOpt1[i,2],4),
                colour = "black",
                vjust = -1
              )
            
            + annotate(
              "point",
              x = round(dfOpt1[i + 1, 2],4) ,
              y = 0,
              colour = "black",
              size = 2
            ) +
              annotate(
                "text",
                x = round(dfOpt1[i + 1, 2],4)+0.15,
                y = 0,
                label =round( dfOpt1[i + 1, 2],4),
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
                x = 0+0.15,
                y = -10,
                label ="Punto mínimo",
                colour = "black",
                vjust = 1.5
              )
          )
        
    }
    
    
  }
  
  dfOpt2 <- newtonII(function(x)
    (3*x^2+8*x),  function(x)
      (6*x+8)
    , 0.75)
  colnames (dfOpt2) <- c ( 'X_n', 'F(X_n)','F``(X_n)', 'X_(n+1)')
  dfOpt1 <- newton(function(x)
    (3*x^2+8*x),  function(x)
      (6*x+8)
    , 0.75)
  xNOpt <- seq(-1.5, 1.5, 0.05)
  
  
  #GRAFICA
  
  xNOptdf <- as.data.frame(xNOpt)
  gOpt <- ggplot(xNOptdf, aes(x = xNOpt, y = (3*xNOpt^2+8*(xNOpt)))) + annotate(
    "text",
    x =1.5,
    y =  16,
    label = "f '(x)",
    colour = "black",
    vjust = 1.5,size=5
  ) +geom_line( aes(x = xNOpt, y =  (xNOpt^3+4*(xNOpt)^2-10)),color="black",size=1)+ annotate(
    "text",
    x =1.5,
    y =  -2,
    label = "f(x)",
    colour = "black",
    vjust = 1.5,size=5
  ) +
    labs(x = "X", y = "F(X)") +
    theme(axis.title.x = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +
    theme(axis.title.y = element_text(
      face = "bold",
      vjust = -0.5,
      colour = "#f2540c",
      size = rel(1)
    )) +theme_bw()+
    geom_line(data = xNOptdf) +
    geom_line(color = "blue", size = 1)
  
  
  #Y=mx+b  m=> SLOPE b=>INTERCEPT
  slopeOpt <- function(x) {
    xincrementO <- 0.2
    yincrementO <-
      (3*(x+0.1)^2+8*(x+0.1))-(3*(x-0.1)^2+8*(x-0.1))
    slopeO <- yincrementO / xincrementO
    return(slopeO)
  }
  
  interceptOpt <- function(x) {
    interceptO <- (3*x^2+8*(x)) - slopeOpt(x) * x
    yvaluesO <- slopeOpt(x) * xNOpt + interceptO
    xNTrigdf <- as.data.frame(cbind(xNOptdf, yvaluesO))
    return(interceptO)
  }
  
  output$graficaOptN <- renderPlot(gOpt)
  
  observeEvent(input$botonOptN, newtonOpt(input$botonOptN, 5))
  
  
  
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
    
      if(tab[,3:tam[2]]== pracma::zeros(tam[1],tam[2])){ 
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
        
          output$IntNew <- DT::renderDataTable({
            v$m
          })
          
      }
    }else{
      polinomio<-polinomioInt(v$m[1,],v$m[,1])
      output$textoIntNew<-renderText(paste("<h4>El polinomio de interpolación es: <br><br> P(x) =", polinomio," <h4>"))
      dataX <- as.data.frame(v$m[,1])
       gInterInf<-ggplot(v$m )+geom_line(aes(x=v$m[,1],y=v$m[,2]))+ theme_bw()+ labs(x = "X", y = "F(X)")
      output$graficaIntNew <-
        renderPlot(gInterInf)
    }
    
   
    
  })
  
  
  
}

shinyApp(ui = ui, server = server)

