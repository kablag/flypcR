library(shiny)
library(ggplot2)

source("R/world.R")
source("R/flywire_data.R")
source("R/flywire_brain.R")
source("R/brain.R")
source("R/controller.R")

ui <- fluidPage(
  titlePanel("FlyPCR — раскапка ПЦР мозгом мухи"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_samples", "Число образцов ДНК", value = 8, min = 0, max = 60, step = 1),
      helpText("Будут добавлены 3 ПКО, ОКО выделения и водный ОКО."),
      actionButton("reset", "Сформировать протокол", class = "btn-primary"),
      actionButton("step", "Один шаг"),
      actionButton("run_all", "Выполнить всё"),
      checkboxInput("play", "Режим просмотра", value = FALSE),
      sliderInput("steps_per_second", "Скорость, шагов в секунду",
                  min = 1, max = 30, value = 8, step = 1),
      hr(),
      verbatimTextOutput("status")
    ),
    mainPanel(
      plotOutput("plate", height = "470px"),
      h4("Операции"),
      tableOutput("operations"),
      h4("Последние действия"),
      verbatimTextOutput("log")
    )
  )
)

server <- function(input, output, session) {
  world <- reactiveVal(new_world(8))
  brain <- reactiveVal(new_brain())

  do_step <- function() {
    w <- world()
    if (w$action_index >= nrow(w$actions)) return(invisible(FALSE))
    action <- w$actions[w$action_index + 1L, ]
    b <- brain_protocol_step(brain(), action)
    world(apply_action(w, brain_action(b)))
    brain(b)
    invisible(TRUE)
  }

  observeEvent(input$step, do_step())

  observe({
    playing <- isTRUE(input$play)
    speed <- input$steps_per_second
    if (!playing || is.null(speed) || speed < 1) return()

    invalidateLater(max(20, round(1000 / speed)), session)
    advanced <- isolate(do_step())
    if (!advanced) updateCheckboxInput(session, "play", value = FALSE)
  })

  observeEvent(input$run_all, {
    while (do_step()) { }
  })

  observeEvent(input$reset, {
    n <- suppressWarnings(as.integer(input$n_samples))
    validate(need(!is.na(n) && n >= 0 && n <= 60,
                  "Введите от 0 до 60 образцов"))
    world(new_world(n))
    brain(new_brain())
    updateCheckboxInput(session, "play", value = FALSE)
  })

  output$plate <- renderPlot({
    w <- world()
    plate <- expand.grid(row = LETTERS[1:8], col = 1:12)
    plate$y <- match(plate$row, LETTERS)
    df <- merge(plate, w$wells, by = c("row", "col"), all.x = TRUE, sort = FALSE)
    df$used <- !is.na(df$reaction)
    df$complete <- !is.na(df$volume) & abs(df$volume - 25) < 1e-8

    df$deck_x <- 10 + df$col - 1
    df$deck_y <- 2 + df$y - 1
    components <- data.frame(
      x = c(1, 3, 5, 7), y = 2,
      label = c("PCR mix\n10X", "Праймеры\n5X", "Вода", "Мастер-\nмикс"), type = "Реагенты")
    equipment <- data.frame(
      x = c(1,3,5,7,10,12,14,16,19), y = 15,
      label = c("P1000","P200","P20","P10","Tips\n1000","Tips\n200","Tips\n20","Tips\n10","Отходы"), type = "Оснащение")
    samples <- data.frame(
      x = 1 + ((seq_len(w$n_reactions) - 1) %% 7),
      y = 5 + ((seq_len(w$n_reactions) - 1) %/% 7),
      label = w$wells$sample, type = "Матрицы")
    objects <- rbind(components, equipment, samples)

    ggplot() +
      geom_point(data = df, aes(deck_x, deck_y, fill = complete, alpha = used),
                 size = 7, shape = 21, colour = "grey35") +
      geom_text(data = subset(df, used), aes(deck_x, deck_y, label = ifelse(volume > 0, volume, "")), size = 2.7) +
      geom_point(data = objects, aes(x, y, colour = type), size = 5, shape = 15) +
      geom_text(data = objects, aes(x, y + 0.45, label = label), size = 2.5, lineheight = .85) +
      geom_point(data = data.frame(x = w$fly$x, y = w$fly$y), aes(x, y),
                 size = 7, shape = 8, stroke = 1.4, colour = "#D62728") +
      geom_text(data = data.frame(x = w$fly$x, y = w$fly$y), aes(x, y - .5), label = "муха", colour = "#D62728", size = 3) +
      scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#62C370"), guide = "none") +
      scale_alpha_manual(values = c(`FALSE` = 0.18, `TRUE` = 1), guide = "none") +
      scale_colour_manual(values = c(Реагенты = "#1F77B4", Оснащение = "#555555", Матрицы = "#9467BD"), guide = "none") +
      scale_x_continuous(limits = c(0, 22), breaks = 10:21, labels = 1:12) +
      scale_y_reverse(limits = c(16, 0), breaks = 2:9, labels = LETTERS[1:8]) + coord_fixed() +
      labs(x = NULL, y = NULL, subtitle = "Планшет справа; красная звезда — положение мухи") +
      theme_minimal(base_size = 13)
  })

  output$status <- renderPrint({
    w <- world(); b <- brain()
    cat("Реакций:", w$n_reactions, "(образцы + 5 контролей)\n")
    cat("Мастер-микс:", 23 * w$n_reactions, "мкл\n")
    cat("Готово:", reward(w), "из", w$n_reactions, "реакций\n")
    cat("Действие:", w$action_index, "из", nrow(w$actions), "\n")
    cat("Шагов выполнено:", w$step_count, "\n")
    cat("Просмотр:", if (isTRUE(input$play))
      paste(input$steps_per_second, "шагов/с") else "остановлен", "\n")
    cat("Положение мухи:", w$fly$x, w$fly$y, "\n")
    cat("Пипетка:", ifelse(is.na(w$mounted_pipette), "нет", w$mounted_pipette),
        "| наконечник:", ifelse(w$has_tip, "да", "нет"),
        "| несёт:", w$carried_volume, "мкл\n")
    cat("Наконечники:\n"); print(w$tips)
    cat("Всего:", sum(w$tips), "(оптимум после полного выполнения:", w$n_reactions + 4, ")\n")
    cat("\nFlyWire v", b$version, ", MDN:\n", sep = "")
    print(round(b$mdn_output, 4))
  })

  output$operations <- renderTable({
    w <- world()
    aggregate(cbind(`Объём, мкл` = volume, Ходов = strokes) ~ stage + source + pipette,
              data = w$tasks, FUN = sum)
  }, striped = TRUE, spacing = "xs", digits = 2)

  output$log <- renderText({
    x <- tail(world()$log, 10)
    if (!length(x)) "Протокол ещё не запущен." else paste(x, collapse = "\n")
  })
}

shinyApp(ui, server)
