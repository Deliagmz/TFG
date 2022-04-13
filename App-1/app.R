library (shiny)
library(shinythemes)
Sys.setlocale(category = "LC_ALL", locale = "Spanish")



##FUNCIONES BIPART
feval <- function(f, ...) do.call(f, list(...))


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
                              I = paste0("[", a, ", ", b, "]"),
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
          
          h3("¿Comenzamos? :)"),
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
                                                             h4("Método de bipartición: "),
                                                             br(),
                                                             helpText("Dada una ecuación \\(f(x)=0\\) con \\(f(x)\\) continua."),
                                                             
                                                             column(
                                                               12,
                                                               offset = 0.5,
                                                               helpText("1) Busco un intervalo \\([a , b ]\\) tal que \\(f(a)f(b)<0\\)"),
                                                               
                                                               helpText("2) Repetir hasta criterio de parada:"),
                                                               column(
                                                                 12,
                                                                 offset = 0.5,
                                                                 helpText("Tomo  \\(x=\\frac{a+b}{2}\\)."),
                                                                 helpText(
                                                                   "Evalúo \\(f(x):\\begin{cases}
                         f(x)=0  & \\text{Parar,} \\\\
                         f(a)f(x)<0 & \\text{Nuevo intervalo [a,x]}\\\\
                         f(x)f(b)<0 & \\text{Nuevo intervalo [x,b]}
                                 \\end{cases}\\! \\) "
                                                                 )
                                                               ),
                                                               helpText("3) Devolver el valor \\(x=\\frac{a+b}{2}\\) como solución. ")
                                                             ),
                                                             
                                                             
                                                           ),
                                                           
                                                           column(
                                                             6,
                                                             offset = 0.5,
                                                             br(),
                                                             br(),
                                                             h4("Método de Newton: "),
                                                             br(),
                                                             helpText("Partimos de un punto inicial \\(x_0\\) e iteramos mediante: "),
                                                             column(
                                                               12,
                                                               offset = 2,
                                                               helpText("\\(x_k+1= x_k -\\frac{f(x_k)}{f'(x_k)} \\)"),
                                                             ),
                                                             
                                                             
                                                           ),
                                                         ))),
        tabPanel(
          "Ejemplos Bipartición",
          uiOutput("Biparticion"),
          mainPanel(
            navlistPanel(
              
              tabPanel("Ecuación exponencial",  mainPanel(
                fluidRow(
                br(),
                h4("Ecuación exponencial con intervalo [0,1]: $$e^{(3x)}-4$$ "),
                br(),
               
                br(),br(), br()
              ),
              fluidRow(splitLayout( style = "height: 200px;
                                          position: relative; ",
                                    style = "text-align: center",
                                   cellWidths = 400,
                                  style = "padding: 40px",
                
                       tableOutput("exp"),
                       plotOutput("graficaExp")),
                
              ),
              fluidRow( actionButton("botonexp", ">>>"),)
              )),
              
              
              
              tabPanel("Ecuación trigonométrica", uiOutput("trig"),
                       
                       actionButton("botontrig",">>>")),
              
              tabPanel("Punto de corte entre funciones", "corte"),
            
              tabPanel("Optimización", uiOutput("opt"),
                       
                       actionButton("botonOpt",">>>")),
              
            )
          ),
          
        ),
        tabPanel(
          "Ejemplos Newton",
          uiOutput("Newton"),
          mainPanel(
            navlistPanel(
              
              tabPanel("Ecuación exponencial", "exp"),
              
              tabPanel("Ecuación trigonométrica", "trig"),
              
              tabPanel("Punto de corte entre funciones", "corte"),
              
              tabPanel("Optimización", "opt"),
              
            )
          ),
        )
      )
    ),
    
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
                             helpText(
                               "Consiste en hallar una función polinómica del menor grado posible que pase por los puntos dados."
                             ),
                             
                             column(
                               12,
                               offset = 0.5,
                               helpText(
                                 "Si consideramos la base {\\(1,x,x^2,...,x^n\\)}  el polinomio \\(P_n (x)\\) que buscamos es
                                                 $$P_n (x)= a_0+a_1 x+a_2 x^2+...+a_n x^n$$"
                               ),
                               
                               helpText(
                                 "Como queremos que \\(P_n(x_i) = f(x_i)\\) con \\(i =0,1,...,n\\), para encontrar los coeficientes
                                                    \\(a_0,a_1,a_2,...,a_n\\) hemos de resolver el sistema de ecuaciones lineales :"
                               ),
                               column(
                                 12,
                                 offset = 0.5,
                                 helpText(
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
                               helpText(
                                 "La matriz de coeficientes anterior es una matriz de Vandermonde, cuyo determinante viene dado por \\( \\prod\\nolimits_{j<k}(x_j-x_k)\\) .
                                           Por tanto, si todos los nodos son distintos, dicho determinante es distinto de
                                            cero, y el sistema anterior tiene una única solución. "
                               ),
                               br(),
                               
                               helpText(
                                 "El error cometido al aproximar\\(f(x)\\) por \\(P_n(x)\\) se puede expresar como  "
                               ),
                               
                               helpText(
                                 "$$ e_n(x) = \\frac{f^{n+1}(\\xi)}{(n+1)!}\\ (x-x_0)(x-x_1)...(x-x_n) $$ "
                               ),
                               helpText("con \\(\\xi\\) entre el mayor y el menor de los nodos. ")
                             ),
                           ),
                         )),),
      
      tabPanel("De Newton", uiOutput("INewton"),
               tabsetPanel(tabPanel(
                 "Teoría", uiOutput("TINewton"),
                 mainPanel(fluidRow
                           (
                             column(
                               12,
                               br(),
                               h4("Interpolación de Newton en diferencias divididas"),
                               br(),
                               helpText("Se llama diferencias divididas a los coeficientes ."),
                               helpText("$$f [ x_j,x_j+1,...,x_i-1,x_i]$$"),
                               helpText("Se pueden calcular usando el esquema:"),
                               column(
                                 12,
                                 offset = 0.5,
                                 helpText(
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
                               strong("Fórmula de Newton"),
                               helpText(
                                 "El polinomio de interpolación se pude escribir como: $$ P_n(x) = f[x_0]+ f[x_0,x_1](x-x_0)+f[x_0,x_1,x_2](x-x_0)(x-x_1)+...+f[x_0,x_1,x_2,...,x_n](x-x_0)(x-x_1)...(x-x_n).$$"
                               ),
                               strong("Error de interpolación"),
                               helpText(
                                 "Si \\(f\\) es una función continua en \\([a , b ], K \\) veces derivables en \\((a , b )\\) y la derivada k-ésima es continua, entonces:
                                 $$ E(x) = \\frac{f^{n+1}(\\xi)}{(n+1)!}\\prod\\limits_{i=0}^n (x-x_i) $$ "
                                 
                               ),
                             )
                           ))
               ),
               tabPanel(
                 "Ejemplos", uiOutput("EINewton"), mainPanel(fluidRow
                                                            (
                                                              column(
                                                                12,
                                                                br(),
                                                                h4("Interpolación de Newton en diferencias divididas"),
                                                                br(),
                                                                numericInput("X", "X:", 0, min = 0, max = 100),
                                                                numericInput("Y", "Y:", 0, min = 0, max = 100),
                                                                plotOutput("grafica", width = "100%", height = "400px", inline = FALSE),
                                                                column(
                                                                  12,
                                                                  offset = 0.5,
                                                                  helpText(
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
                                                                strong("Fórmula de Newton"),
                                                                helpText(
                                                                  "El polinomio de interpolación se pude escribir como: $$ P_n(x) = f[x_0]+ f[x_0,x_1](x-x_0)+f[x_0,x_1,x_2](x-x_0)(x-x_1)+...+f[x_0,x_1,x_2,...,x_n](x-x_0)(x-x_1)...(x-x_n).$$"
                                                                ),
                                                                strong("Error de interpolación"),
                                                                helpText(
                                                                  "Si \\(f\\) es una función continua en \\([a , b ], K \\) veces derivables en \\((a , b )\\) y la derivada k-ésima es continua, entonces:
                                 $$ E(x) = \\frac{f^{n+1}(\\xi)}{(n+1)!}\\prod\\limits_{i=0}^n (x-x_i) $$ "
                                                                  
                                                                ),
                                                              )
                                                            ))),
              
              
                 )), 
      tabPanel("Inversa", uiOutput("IInversa"), 
               tabsetPanel(tabPanel(
        "Teoría", uiOutput("TIInversa"),
        mainPanel(fluidRow
                  (
                    column(
                      12,
                      br(),
                      h4("Interpolación inversa"),
                      br(),
                      
                      helpText("Dados unos valores de la función \\(\\left\\{(x_i,y_i)\\right\\}\\) y un valor \\(y\\), el objetivo ahora es encontrar el valor de \\(x\\), tal que \\(f(x)=y\\)  ."),
                      column(
                        12,
                        offset = 0.5,br(),
                        strong("Método 1:"),
                        helpText(
                          "Si los valores \\(y_i\\) son diferentes entre sí, podemos cambiar los papeles de \\(x_i\\) e \\(y_i\\). Lo que hacemos es interpolar a los datos \\((x_i,y_i)\\). "
                        ),br(),
                        strong("Método 2:"),
                        helpText(
                          "Interpolamos una función a los datos \\(\\left\\{(x_i,y_i)\\right\\}\\), de manera que se obtiene una estimación de \\(f(x)\\).Después resolvemos la ecuación \\(f(x)=y\\)  . "
                        )
                      ),
                    )
                  ))
      ),
      tabPanel(
        "Ejemplos", uiOutput("EIInversa")),
      ),
      ),
      tabPanel("Osculatoria",
               uiOutput("IOculatoria"),
               tabsetPanel(tabPanel(
                 "Teoría", uiOutput("TIOsculatoria"),
                 mainPanel(fluidRow
                           (
                             column(
                               12,
                               br(),
                               h4("Interpolación osculatoria"),
                               br(),
                               
                               helpText(
                                 "Este método se utiliza si además de exigir que el polinomio interpolador coincida con la función en los nodos, se impone que ciertas derivadas en dichos nodos también coincidan con las 
                                 correspondientes deirvadas de la función."
                               ),
                               column(
                                 12,
                                 offset = 0.5,
                                 
                                 helpText(
                                   "Utilizaremos el esquema en diferencias divididas visto en la interpolación de Newton, pero repitiendo \\(k_i\\) veces cada nodo \\(x_i\\) y utilizando que:
                                   $$ f[x_i,x_i]=f'(x_i), f[x_i,x_i,x_i]=\\frac{f''(x_i)}{2},...,f[x_i,x_i,...^{n+1},x_i]=\\frac{f^{n)}(x_i)}{n!} $$" ),
                                 br(),
                                
                               ),
                               strong("Interpolacion de Hermite"),
                               helpText("Hablamos de este caso particular cuando se usa como informaión el valor de la función y de su derivada en cada uno de los nodos.
                                        Se trata del polinomio de interpolación osculatoria con \\(k_i=1, i = 0,...,n\\), y de grado por tanto \\(2n+1\\). "),
                               helpText("La fórmula de error en este caso es: $$ e_{2n+1}(x) = \\frac{f^{2n+2}(\\xi)}{(2n+2)!}\\ (x-x_0)^2(x-x_1)^2...(x-x_n)^2$$")
                             )
                           ))
               ),
               tabPanel(
                 "Ejemplos", uiOutput("EIOsculatoria")),
               )
               )
    ),
    tabPanel("Derivación", uiOutput("opcion3")),
    tabPanel("Integración", uiOutput("opcion4")),
    
  )
)


server <- function(input, output) { 
 
  ##BIPARTICION
      #EXPONENCIAL
  bipartExp<-function(i,n,ab){
    if(i<=n){
      output$exp<-renderTable(df[1:i,],bordered = TRUE,colnames=TRUE )
      output$graficaExp<-renderPlot(plot(x,exp(3*x)-4,ylim=c(-5,10),xlim = c(ab[i,1],ab[i,2]),type="l",col="orange",lwd=3))
    }
   
  }
  
  df<-bipart_table(function(x) (exp(3*x)-4),0,1,5,1)
  ab<-bipart(function(x) (exp(3*x)-4),0,1,5)
  
  output$graficaExp<-renderPlot(plot(x,exp(3*x)-4,ylim=c(-5,10),xlim = c(0,1),type="l",col="orange",lwd=3))
  
  observeEvent(input$botonexp,bipartExp(input$botonexp,5,ab))
  
      #TRIGONOMETRICA
  bipartTrig<-function(i,n){
    if(i<=n){
      output$trig<-renderTable(dfT[1:i,])
    }
  }
  dfT<-bipart_table(function(x) sin(x)+cos(x^2),-1,1,5,1)
  
  #output$trig<-renderTable(dfT[1:1,])
  observeEvent(input$botontrig,bipartTrig(input$botontrig,5))
  
    #OPTIMIZACION
  bipartOpt<-function(i,n){
    if(i<=n){
      output$opt<-renderTable(dfO[1:i,])
    }
  }
  dfO = bipart_table(function(x) (x^3-x-1),1,2,6,1)
  output$opt<-renderTable(dfO[1:1,])
  observeEvent(input$botonOpt,bipartOpt(input$botonOpt+1,6))
  
  output$grafica <- renderPlot({
    plot(mtcars$wt, mtcars$mpg)
  })
}

shinyApp(ui = ui, server = server)



