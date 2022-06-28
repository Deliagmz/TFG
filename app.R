library (shiny)
library(shinythemes)
library(ggplot2)
Sys.setlocale(category = "LC_ALL", locale = "Spanish")

##FUNCIONES BIPART
feval <- function(f, ...) do.call(f, list(...))

intNewton_table <- function(x, y) {
  
  n <- length(x)
  print(n)
  z <- pracma::zeros(n,  n+1)
  
  z[, 1] <- t(x)
  z[, 2] <- t(y)
  print(z)
  for (col in 3:(n+1)) {
    
    for (fil in 1:(n-1)) {
      
      z[fil, col] <- (z[fil + 1, col - 1] - z[fil,
                                              col - 1])/(z[fil + col - 2] - z[fil]) 
    }
    
  }
  return(z)
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
                  cellArgs = list(style = "border: 1px solid #f2540c;padding:0.5em;"),
                  cellWidths = 400,
                  tableOutput("exp"),
                  plotOutput("graficaExp")
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botonexp", ">>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
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
                  cellArgs = list(style = "border: 1px solid #f2540c;padding:0.5em;"),
                  cellWidths = 400,
                  tableOutput("trig"),
                  plotOutput("graficaTrig")
                ),
                br(),
                column(
                  width = 12,
                  offset = 8,
                  align = "right",
                  actionButton("botontrig", ">>>", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
                ),
                
              ),
              
            )),
            
            tabPanel("Punto de corte entre funciones", "corte"),
            
            tabPanel("Optimización", mainPanel(
              fluidRow(
                br(),
                h4(
                  "Ecuación de optimización con intervalo [0,2]: $$f(x)= (x^3-x-1)$$ $$f'(x)= (3x^2-1)$$"
                ),
                br(),
              ),
              fluidRow(
                splitLayout(
                  cellArgs = list(style = "border: 1px solid #f2540c;padding:0.5em;"),
                  cellWidths = 400,
                  tableOutput("opt"),
                  plotOutput("graficaOpt")
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
                h4("Ecuación exponencial  $$f(x)=e^{(3x)}+x^2-4$$ "),
                br(),
              ),
              fluidRow(
                splitLayout(
                  cellArgs = list(style = "border: 1px solid #f2540c;padding:0.5em;"),
                  cellWidths = 400,
                  tableOutput("expN"),
                  plotOutput("graficaExpN")
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
            
            tabPanel("Ecuación trigonométrica", "trig"),
            
            tabPanel("Punto de corte entre funciones", "corte"),
            
            tabPanel("Optimización", "opt"),
            
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
                cellArgs = list(style = "border: 1px solid #f2540c;padding:0.5em;"),
                cellWidths = 700,
                DT::dataTableOutput("IntNew"),
                #tableOutput("IntNew"),
                plotOutput("graficaIntNew")
              ),
              br(),
              column(
                width = 12,
                offset = 7.5,
                align = "centre",
                p("Introduce los puntos 2 a 2 : "),
                numericInput("x", "X:", 0, min = -50, max = 50),
                numericInput("y", "Y:", 0, min = -50, max = 50),
                actionButton("botonIntNewP", "Introducir punto"),
                actionButton("botonIntNew", "Sig. Paso", style = "color: #fff; background-color: #f2540c; border-color: #f2540c")
              ),
              
            ),
            
          )),
        )
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
                                   "Hablamos de este caso particular cuando se usa como informaión el valor de la función y de su derivada en cada uno de los nodos.
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
          tabPanel("Teoría", uiOutput("teoría"),
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
          tabPanel("Teoría", uiOutput("teoría"),
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
          tabPanel("Teoría", uiOutput("teoría"),
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
            uiOutput("EDerivacion"),
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

v <- reactiveValues()
v$df <- data.frame(x = numeric(), y = numeric())

server <- function(input, output) { 
  
  ##BIPARTICION
  #EXPONENCIAL
  bipartExp <- function(i, n, ab) {
    if (i <= n) {
      output$exp <- renderTable(df[1:i, ])
      output$graficaExp <-
        renderPlot(plot(
          x,
          exp(3 * x) - 4,
          ylim = c(-5, 10),
          xlim = c(ab[i, 1], ab[i, 2]),
          type = "l",
          col = "#f2540c",
          lwd = 3
        ))
    }
    
  }
  x<-seq(0,1,0.05)
  df <- bipart_table(function(x)
    (exp(3 * x) - 4), 0, 1, 5, 1)
  ab <- bipart(function(x)
    (exp(3 * x) - 4), 0, 1, 5)
  
  output$graficaExp <-
    renderPlot(plot(
      x,
      exp(3 * x) - 4,
      ylim = c(-5, 10),
      xlim = c(0, 1),
      type = "l",
      col = "#f2540c",
      lwd = 3
    ))
  
  observeEvent(input$botonexp, bipartExp(input$botonexp, 5, ab))
  
  #TRIGONOMETRICA
  bipartTrig <- function(i, n, abT) {
    if (i <= n) {
      output$trig <- renderTable(dfT[1:i, ])
      output$graficaTrig <-
        renderPlot(plot(
          xT,
          sin(xT) + cos(xT ^ 2),
          ylim = c(-5, 10),
          xlim = c(abT[i, 1], abT[i, 2]),
          type = "l",
          col = "#f2540c",
          lwd = 3
        ))
    }
  }
  xT<-seq(-1,1,0.05)
  dfT <- bipart_table(function(xT)
    sin(xT) + cos(xT ^ 2), -1, 1, 5, 1)
  abT <- bipart(function(xT)
    sin(xT) + cos(xT ^ 2), -1, 1, 5)
  
  output$graficaTrig <-
    renderPlot(plot(
      xT,
      sin(xT) + cos(xT ^ 2),
      ylim = c(-5, 10),
      xlim = c(-1, 1),
      type = "l",
      col = "#f2540c",
      lwd = 3
    ))
  
  observeEvent(input$botontrig, bipartTrig(input$botontrig, 5, abT))
  
  #OPTIMIZACION
  bipartOpt <- function(i, n, abO) {
    if (i <= n) {
      output$opt <- renderTable(dfO[1:i, ])
      output$graficaOpt <-
        renderPlot(plot(
          xO,
          (3 * xO ^ 2 - 1),
          ylim = c(-5, 10),
          xlim = c(abO[i, 1], abO[i, 2]),
          type = "l",
          col = "#f2540c",
          lwd = 3
        ))
    }
  }
  xO<-seq(0,2,0.05)
  dfO <- bipart_table(function(xO)
    (3 * xO ^ 2 - 1), 0, 2, 6, 1)
  abO <- bipart(function(xO)
    (3 * xO ^ 2 - 1), 0, 2, 6)
  
  output$graficaOpt <-
    renderPlot(plot(
      xO,
      (3 * xO ^ 2 - 1),
      ylim = c(-5, 10),
      xlim = c(0, 2),
      type = "l",
      col = "#f2540c",
      lwd = 3
    ))
  
  observeEvent(input$botonopt, bipartOpt(input$botonopt, 6, abO))
  
  ##NEWTON
  #EXPONENCIAL
  newtonExp <- function(i, n) {
    if (i <= n) {
      output$expN <- renderTable(dfN[1:i,])
      output$graficaExpN <-
        renderPlot(
          g +  geom_abline(
            slope = slope(dfN[i, 2]),
            intercept = intercept(dfN[i, 2]),
            color = "#f2540c",
            size = 1
          )
          + annotate(
            "point",
            x = dfN[i, 2] ,
            y = dfN[i, 4],
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = dfN[i, 2],
              y = dfN[i, 4],
              label = dfN[i, 2],
              colour = "black",
              vjust = -1
            )
          
          + annotate(
            "point",
            x = dfN[i + 1, 2] ,
            y = 0,
            colour = "black",
            size = 2
          ) +
            annotate(
              "text",
              x = dfN[i + 1, 2],
              y = 0,
              label = dfN[i + 1, 2],
              colour = "black",
              vjust = 1.5
            )
        )
      
    }
    
    
  }
  
  dfN <- newton(function(x)
    (exp(3 * x) + x ^ 2 - 4),  function(x)
      (3 * exp(3 * x) + 2 * x)
    , 4)
  xN <- seq(0, 4, 0.05)
  
  
  #GRAFICA
  
  xNT <- as.data.frame(xN)
  g <- ggplot(xNT, aes(x = xN, y =  exp(3 * xN) + xN ^ 2 - 4)) +
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
    )) +
    geom_line(data = xNT) +
    geom_line(color = "blue", size = 1)
  
  
  #Y=mx+b  m=> SLOPE b=>INTERCEPT
  slope <- function(x) {
    xincrement <- 0.2
    yincrement <-
      exp(3 * (x + 0.1)) + (x + 0.1) ^ 2 - 4 - exp(3 * (x - 0.1)) + (x - 0.1) ^
      2 - 4
    slope <- yincrement / xincrement
    return(slope)
  }
  
  intercept <- function(x) {
    intercept <- exp(3 * x) + x ^ 2 - 4 - slope(x) * x
    yvalues <- slope(x) * xN + intercept
    xNT <- as.data.frame(cbind(xNT, yvalues))
    return(intercept)
  }
  
  output$graficaExpN <- renderPlot(g)
  
  observeEvent(input$botonexpN, newtonExp(input$botonexpN, 10))
  
  
  ##INTERPOLACION NEWTON
  
  
  
  ##INTRODUCIR DATOS
  observeEvent(input$botonIntNewP,{
    req(input$x,input$y) 
    
    
    tmp <- data.frame(x = input$x,y = input$y)
    
    ##HACER PLOT REACTIVE## 
    
    #   g <- ggplot(data = tmp) + 
    #   geom_point(aes(tmp[,1], tmp[,2]))+
    #     annotate(
    #   "point",
    #   x = input$x,
    #   y = input$y,
    #   colour = "blue",
    #   size = 2
    # ) +
    #   annotate(
    #     "text",
    #     x = input$x,
    #     y = input$y,
    #     label = input$x,
    #     colour = "black",
    #     vjust = -1
    #   )
    #   
    v$df <- rbind(v$df,tmp)
    output$IntNew <- DT::renderDataTable({
      v$df
    }) 
    g<-ggplot((v$df)) +
      geom_point(aes(x =v$df[,1],y=v$df[,2])) 
    output$graficaIntNew <-
      renderPlot(g)
  })
  
  ##REALIZAR TABLA
  
  ##ELIMINAR DATOS QUE SOBRAN
  observeEvent(input$botonIntNew,{
    tabla<- as.data.frame(v$df)
    tab<- intNewton_table(tabla[,1],tabla[,2])
    fx<- tab[,input$botonIntNew+2]
    v$df <- cbind(v$df,
                  fx)
    output$IntNew <- DT::renderDataTable({
      v$df
    })
    
  })
  
  
  
}

shinyApp(ui = ui, server = server)