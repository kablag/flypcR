library(shiny)
library(ggplot2)

source("R/world.R")
source("R/flywire_data.R")
source("R/flywire_brain.R")
source("R/brain.R")
source("R/controller.R")

ui <- fluidPage(

  titlePanel("FlyPCR"),

  fluidRow(

    column(
      8,
      plotOutput("plate")
    ),

    column(
      4,

      h3("Fly brain"),

      verbatimTextOutput("status"),

      actionButton(
        "step",
        "Brain step"
      ),

      actionButton(
        "reset",
        "Reset"
      )
    )
  )
)


server <- function(input, output, session) {

  world <- reactiveVal(new_world())
  brain <- reactiveVal(new_brain())

  observeEvent(input$step, {

    w <- world()
    b <- brain()

    obs <- observe_world(w)

    b <- brain_step(
      b,
      obs
    )

    action <- brain_action(b)

    w <- apply_action(
      w,
      action
    )

    brain(b)
    world(w)
  })


  observeEvent(input$reset, {

    world(new_world())
    brain(new_brain())

  })


  output$plate <- renderPlot({

    w <- world()

    df <- w$wells

    df$y <- match(
      df$row,
      LETTERS
    )

    ggplot(
      df,
      aes(col, y)
    ) +

      geom_point(
        aes(size = volume),
        shape = 21
      ) +

      geom_point(
        data = data.frame(
          col = w$pipette$x,
          y = w$pipette$y
        ),
        size = 8,
        shape = 4
      ) +

      scale_y_continuous(
        breaks = 1:8,
        labels = LETTERS[1:8]
      ) +

      scale_x_continuous(
        breaks = 1:12
      ) +

      coord_fixed() +

      labs(
        x = NULL,
        y = NULL
      ) +

      theme_minimal()

  })


  output$status <- renderPrint({

    w <- world()
    b <- brain()

    obs <- observe_world(w)

    cat(
      "Pipette:",
      w$pipette$x,
      w$pipette$y,
      "\n"
    )

    cat(
      "Sensors:\n"
    )

    print(obs)

    cat(
      "\nMotor output:\n"
    )

    print(
      brain_action(b)
    )

    cat("\nFlyWire v", b$version, " MDN activity:\n", sep = "")
    print(round(b$mdn_output, 4))

    cat(
      "\nReward:",
      reward(w),
      "\n"
    )
  })
}


shinyApp(ui, server)
