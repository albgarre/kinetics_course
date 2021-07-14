
library(shiny)
library(shinyWidgets)
library(bs4Dash)

library(tidyverse)
library(biogrowth)

library(thematic)
thematic_shiny()

## UI---------------------------------------------------------------------------

ui <- dashboardPage(
  dashboardHeader(title = "Variability, uncertainty and Monte Carlo"),
  footer = dashboardFooter(
    left = a(
      href = "https://github.com/albgarre/kinetics_course",
      target = "_blank", "GitHub page"
    ),
    right = "Reaction kinetics in food science"
  ),
  dashboardSidebar(disable=TRUE),
  dashboardBody(
    fluidRow(
      column(4,
             bs4Card(title = "Input parameters", width = 12,
                     footer = actionBttn("go", "Simulate!"),
                     fluidRow(
                       column(12, h4("delta (min)"))
                     ),
                     fluidRow(
                       column(6, numericInput("mean_delta", "Mean", 1, min = 0)),
                       column(6, numericInput("sd_delta", "SD", .1, min = 0))
                     ),
                     fluidRow(
                       column(12, h4("beta (·)"))
                     ),
                     fluidRow(
                       column(6, numericInput("mean_beta", "Mean", 1.25, min = 0)),
                       column(6, numericInput("sd_beta", "SD", .2, min = 0))
                     ),
                     fluidRow(
                       column(12, h4("log N0 (log CFU/g)"))
                     ),
                     fluidRow(
                       column(6, numericInput("mean_logN0", "Mean", 8, min = 0)),
                       column(6, numericInput("sd_logN0", "SD", .5, min = 0))
                     ),
                     fluidRow(column(12, numericInput("max_time", "Treatment time (min)", 3, min = 0))),
                     fluidRow(column(12, numericInput("niter", "Number of simulations", 10, min = 0))),
                     hr(),
                     fluidRow(
                       column(6, numericInput("seed", "Seed", 12421, min = 0)),
                       column(6, actionBttn("reset_seed", "Reset seed"))
                     )
             )
      ),
      column(8,
             fluidRow(
               bs4Card(
                 title = "Input distribution",
                 plotOutput("par_distribution")
               ),
               bs4Card(
                 title = "Model predictions",
                 plotOutput("plot_curves")
               ),
             ),
             fluidRow(
               bs4Card(
                 title = "Distribution of the survivors",
                 plotOutput("surv_distribution")
               ),
               bs4Card(
                 title = "Numerical results",
                 tableOutput("surv_table")
               )
             ),
             fluidRow(
               bs4Card(
                 title = "Details",
                 collapsed = TRUE,
                 tableOutput("calculations")
               )
             )
      )
    )

  )
)

## server ----------------------------------------------------------------------
server <- function(input, output) {

  ## Reactive values

  micropars <- reactiveVal()
  my_sims <- reactiveVal()

  ## Resetting the seed

  observeEvent(input$reset_seed, {
    print(paste("Seed set to", input$seed))
    set.seed(input$seed)
  })

  ## Calculations

  observeEvent(input$go, withProgress(message = "Simulating...", {

    ## Growth curves

    pars <- tibble(
      logN0 = rnorm(input$niter, input$mean_logN0, input$sd_logN0),
      delta = rnorm(input$niter, input$mean_delta, input$sd_delta),
      beta = rnorm(input$niter, input$mean_beta, input$sd_beta)
    ) %>%
      mutate(niter = row_number())

    micropars(pars)


    out <- pars %>%
      split(.$niter) %>%
      map(.,
          ~ tibble(time = seq(0, input$max_time, length = 1000),
                   logN = .$logN0 - (time/.$delta)^.$beta
          )
      )

    my_sims(out)

  }))

  ## Plot of the inactivation curves

  output$plot_curves <- renderPlot({

    validate(
      need(my_sims(), "")
    )

    my_sims() %>%
      imap_dfr(., ~ mutate(.x, sim = .y)) %>%
      ggplot() +
      geom_line(aes(x = time, y = logN, colour = sim)) +
      theme(legend.position = "none")

  })

  ## Plot of the input parameters

  output$par_distribution <- renderPlot({

    validate(need(micropars(), ""))

    micropars() %>%
      select(-niter) %>%
      pivot_longer(everything(), names_to = "par", values_to = "value") %>%
      ggplot() +
      geom_histogram(aes(value)) +
      facet_wrap(~par, scales = "free")

  })

  ## Plot of the survivors

  output$surv_distribution <- renderPlot({

    validate(
      need(my_sims(), "")
    )

    my_sims() %>%
      map_dfr(., ~ tail(., 1)) %>%
      ggplot() +
      geom_histogram(aes(logN))

  })

  ## Table of indexes

  output$surv_table <- renderTable({

    validate(
      need(my_sims(), "")
    )

    my_sims() %>%
      map_dfr(., ~ tail(., 1)) %>%
      summarize(`Mean` = mean(logN, na.rm = TRUE),
                `Median` = median(logN, na.rm = TRUE),
                Q90 = quantile(logN, .9, na.rm = TRUE),
                Q99 = quantile(logN, .99, na.rm = TRUE),
                `Q99.9` = quantile(logN, .999, na.rm = TRUE))

  })

  output$calculations <- renderTable({

    validate(
      need(micropars(), "")
    )

    micropars() %>%
      mutate(logNf = logN0 - (input$max_time/delta)^beta) %>%
      select(niter, everything())

  })


}

# Run the application
shinyApp(ui = ui, server = server)
