# Universidad del Valle - Especialización en Estadística Aplicada
# Proyecto Final - Técnicas de Muestreo
# Sebastián Tutistar Valencia - Miguel Angel

library(tidyverse)
library(janitor)
library(lubridate)
library(skimr)

# Cargar datos
df <- SECOP_Contratos_Adj_2020.2024 %>% clean_names()

# Selección de columnas
df_sel <- df %>%
  select(
    id_del_proceso, entidad, nombre_del_proveedor_adjudicado,
    descripcion_del_procedimiento, modalidad_de_contratacion,
    proveedores_invitados, proveedores_con_invitacion_directa,
    proveedores_unicos_con_respuestas, valor_total_adjudicacion,
    departamento_entidad, ciudad_entidad, adjudicado,
    estado_del_procedimiento, estado_resumen,
    fecha_de_publicacion_del_proceso, fecha_adjudicacion, duracion
  )

# Normalización
df_sel <- df_sel %>%
  mutate(
    fecha_de_publicacion_del_proceso = dmy(trimws(fecha_de_publicacion_del_proceso)),
    fecha_adjudicacion = dmy(trimws(fecha_adjudicacion)),
    valor_total_adjudicacion = as.numeric(gsub("[$,]", "", valor_total_adjudicacion)),
    proveedores_invitados = as.numeric(proveedores_invitados),
    duracion = as.numeric(duracion)
  )

# Población
df_poblacion <- df_sel %>% filter(!is.na(departamento_entidad), !is.na(modalidad_de_contratacion), departamento_entidad != "No Definido")
N <- nrow(df_poblacion)
H <- n_distinct(df_poblacion$departamento_entidad)

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("AVANCE 2: DISEÑO Y TAMAÑO DE MUESTRA\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# 1. CÁLCULO TAMAÑO DE MUESTRA
e <- 0.05; alfa <- 0.05; z <- 1.96; p <- 0.3; DEFF <- 0.85; nr <- 0.10
n0 <- (z^2 * p * (1-p)) / (e * p)^2
n_final <- ceiling((n0 * DEFF * N) / (n0 * DEFF + N - 1) / (1 - nr))

cat("Tamaño de muestra:\n")
cat("• e=5%, confianza=95%, DEFF=0.85, NR=10%\n")
cat("• n =", n_final, "\n\n")

# 2. COMPARACIÓN DE DISEÑOS
set.seed(123)

# MAS
muestra_mas <- df_poblacion %>% slice_sample(n = n_final)
dom_mas <- muestra_mas %>% count(departamento_entidad, modalidad_de_contratacion) %>% filter(n>=5) %>% nrow()

# ESTRATIFICADO
afijacion <- df_poblacion %>%
  count(departamento_entidad, name = "N_h") %>%
  mutate(W_h = N_h/N, n_h = pmax(2, pmin(round(n_final * W_h), N_h)))
afijacion$n_h[which.max(afijacion$N_h)] <- afijacion$n_h[which.max(afijacion$N_h)] + (n_final - sum(afijacion$n_h))

muestra <- df_poblacion %>%
  split(.$departamento_entidad) %>%
  map2_df(afijacion$n_h, ~slice_sample(.x, n = .y)) %>%
  left_join(afijacion, by = "departamento_entidad") %>%
  mutate(peso = N_h/n_h)

dom_est <- muestra %>% count(departamento_entidad, modalidad_de_contratacion) %>% filter(n>=5) %>% nrow()

cat("Comparación:\n")
cat("                    MAS    Estratificado\n")
cat("Cobertura deptos:  ", n_distinct(muestra_mas$departamento_entidad), "/", H, "      ", H, "/", H, "\n")
cat("Dominios (n≥5):    ", dom_mas, "          ", dom_est, "\n")
cat("DEFF:               1.0          0.85\n")
cat("\nElegido: ESTRATIFICADO (mejor cobertura y precisión)\n\n")

cat("═══════════════════════════════════════════════════════════════\n")
cat("AVANCE 3: EJECUCIÓN Y ESTIMACIÓN\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# 3. CÁLCULO IRI
Q1 <- quantile(muestra$valor_total_adjudicacion, 0.25, na.rm=T)
Q3 <- quantile(muestra$valor_total_adjudicacion, 0.75, na.rm=T)
IQR <- Q3 - Q1

muestra <- muestra %>%
  mutate(
    flag1 = as.integer(modalidad_de_contratacion == "Contratación directa"),
    flag2 = as.integer(is.na(proveedores_unicos_con_respuestas) | proveedores_unicos_con_respuestas <= 1),
    flag3 = as.integer(!is.na(proveedores_con_invitacion_directa) & proveedores_con_invitacion_directa > 0),
    flag4 = as.integer(!is.na(duracion) & (duracion < 5 | duracion > 365)),
    flag5 = as.integer(!is.na(valor_total_adjudicacion) & 
                         (valor_total_adjudicacion < (Q1-1.5*IQR) | valor_total_adjudicacion > (Q3+1.5*IQR))),
    IRI = 0.25*flag1 + 0.25*flag2 + 0.15*flag3 + 0.15*flag4 + 0.20*flag5
  )

# 4. ESTIMACIÓN GLOBAL
est_estratos <- muestra %>%
  filter(!is.na(IRI)) %>%
  group_by(departamento_entidad) %>%
  summarise(y_h = mean(IRI), s2_h = var(IRI), n_h = n(), .groups="drop") %>%
  left_join(afijacion %>% select(departamento_entidad, W_h, N_h), by = "departamento_entidad") %>%
  mutate(contribucion = W_h * y_h, var_contrib = W_h^2 * s2_h/n_h * (1 - n_h/N_h))

IRI_global <- sum(est_estratos$contribucion)
SE <- sqrt(sum(est_estratos$var_contrib))
CV <- (SE/IRI_global)*100

cat("Estimación Nacional:\n")
cat("• IRI:  ", round(IRI_global, 4), " (", round(IRI_global*100, 2), "%)\n")
cat("• SE:   ", round(SE, 5), "\n")
cat("• CV:   ", round(CV, 2), "%\n")
cat("• IC95: [", round(IRI_global-1.96*SE, 4), ", ", round(IRI_global+1.96*SE, 4), "]\n\n")

# 5. ESTIMACIÓN POR DOMINIOS - SOLUCIÓN AL PROBLEMA CV=NA
cat("ESTIMACIÓN POR DOMINIOS (Departamento × Modalidad)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Estrategia: Calcular CV usando variabilidad DENTRO del dominio
est_dominios <- muestra %>%
  filter(!is.na(IRI)) %>%
  group_by(departamento_entidad, modalidad_de_contratacion) %>%
  summarise(
    n_d = n(),
    media_IRI = mean(IRI),
    sd_IRI = sd(IRI, na.rm = TRUE),
    min_IRI = min(IRI),
    max_IRI = max(IRI),
    q25_IRI = quantile(IRI, 0.25, na.rm = TRUE),
    q75_IRI = quantile(IRI, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_d >= 5) %>%
  mutate(
    # Estrategia 1: Si sd > 0, calcular CV normal
    # Estrategia 2: Si sd = 0 (todos iguales), usar rango como proxy de variabilidad
    rango_IRI = max_IRI - min_IRI,
    
    # SE usando desviación estándar tradicional
    SE_tradicional = sd_IRI / sqrt(n_d),
    
    # CV tradicional (puede dar NA si sd=0)
    CV_tradicional = ifelse(media_IRI > 0, (sd_IRI / media_IRI) * 100, NA),
    
    # CV alternativo usando rango cuando sd=0
    # Fórmula: CV_rango ≈ (Rango/4) / Media * 100
    # Justificación: Rango ≈ 4×SD en distribuciones normales
    CV_rango = ifelse(media_IRI > 0, (rango_IRI / 4 / media_IRI) * 100, NA),
    
    # CV final: usar tradicional si existe, sino usar rango
    CV_final = ifelse(!is.na(CV_tradicional), CV_tradicional, CV_rango),
    
    # Método usado
    metodo_CV = case_when(
      !is.na(CV_tradicional) & CV_tradicional > 0 ~ "SD",
      rango_IRI > 0 ~ "Rango",
      TRUE ~ "Constante"
    ),
    
    # Intervalos de confianza
    IC_inf = pmax(0, media_IRI - 1.96 * SE_tradicional),
    IC_sup = pmin(1, media_IRI + 1.96 * SE_tradicional),
    
    # Amplitud del IC como medida alternativa de precisión
    amplitud_IC = IC_sup - IC_inf,
    
    # Precisión categórica basada en CV o amplitud IC
    precision = case_when(
      metodo_CV == "Constante" ~ "★★★★ Perfecta",
      !is.na(CV_final) & CV_final <= 15 ~ "★★★ Excelente",
      !is.na(CV_final) & CV_final <= 30 ~ "★★ Buena",
      !is.na(CV_final) & CV_final <= 50 ~ "★ Aceptable",
      amplitud_IC <= 0.2 ~ "★★ Buena (IC)",
      TRUE ~ "Baja"
    )
  ) %>%
  arrange(desc(media_IRI))

cat("Metodología para manejo de CV:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("• Si SD > 0:     CV tradicional (SD/Media × 100)\n")
cat("• Si SD = 0:     CV basado en rango (Rango/4/Media × 100)\n")
cat("• Si Rango = 0:  IRI constante → Precisión perfecta\n")
cat("• Alternativa:   Amplitud de IC como medida de precisión\n\n")

cat("Dominios estimables (n≥5):", nrow(est_dominios), "\n\n")

cat("Distribución de métodos de CV:\n")
table_metodos <- table(est_dominios$metodo_CV)
for(i in seq_along(table_metodos)) {
  cat("• ", names(table_metodos)[i], ": ", table_metodos[i], " dominios\n", sep="")
}
cat("\n")

cat("Distribución de precisión:\n")
table_precision <- table(est_dominios$precision)
for(i in seq_along(table_precision)) {
  cat("• ", names(table_precision)[i], ": ", table_precision[i], " dominios\n", sep="")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("TOP 15 DOMINIOS CON MAYOR RIESGO\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

top_dominios <- est_dominios %>% 
  head(15) %>%
  transmute(
    `#` = row_number(),
    Depto = str_trunc(departamento_entidad, 16),
    Modal = str_trunc(modalidad_de_contratacion, 22),
    n = n_d, 
    IRI = round(media_IRI, 3),
    SD = round(sd_IRI, 3),
    Rango = round(rango_IRI, 3),
    `CV%` = round(CV_final, 1),
    Método = metodo_CV,
    `AmplIC` = round(amplitud_IC, 3),
    Prec = precision
  )

print(top_dominios, n = 15)

cat("\n\nINTERPRETACIÓN DE RESULTADOS:\n")
cat("─────────────────────────────────────────────────────────────────\n")
cat("• CV% = Coeficiente de variación (menor es mejor)\n")
cat("• Método = Cómo se calculó el CV (SD o Rango)\n")
cat("• AmplIC = Amplitud del IC 95% (menor es más preciso)\n")
cat("• Cuando SD=0 y Rango=0: Todos los contratos tienen el mismo IRI\n")
cat("  → Precisión PERFECTA (no hay variabilidad)\n\n")

# 6. ANÁLISIS DETALLADO DE DOMINIOS CON SD=0
dominios_constantes <- est_dominios %>% 
  filter(metodo_CV == "Constante") %>%
  select(departamento_entidad, modalidad_de_contratacion, n_d, media_IRI)

if(nrow(dominios_constantes) > 0) {
  cat("DOMINIOS CON IRI CONSTANTE (SD=0):\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Estos dominios tienen todos los contratos con el mismo IRI:\n\n")
  print(dominios_constantes %>%
          transmute(
            Depto = str_trunc(departamento_entidad, 20),
            Modal = str_trunc(modalidad_de_contratacion, 25),
            n = n_d,
            IRI = round(media_IRI, 3)
          ), n = Inf)
  cat("\nInterpretación: Homogeneidad total → Máxima confiabilidad en la estimación\n\n")
}

# 7. AJUSTE NO RESPUESTA
set.seed(456)
muestra_resp <- muestra %>% mutate(resp = rbinom(n(), 1, 0.92)) %>% filter(resp == 1)
ajuste <- muestra_resp %>%
  count(departamento_entidad, name = "n_resp") %>%
  left_join(afijacion %>% select(departamento_entidad, n_h), by = "departamento_entidad") %>%
  mutate(factor = n_h/n_resp)

muestra_aj <- muestra_resp %>%
  left_join(ajuste %>% select(departamento_entidad, factor), by = "departamento_entidad") %>%
  mutate(peso_aj = peso * factor)

IRI_aj <- muestra_aj %>%
  filter(!is.na(IRI)) %>%
  group_by(departamento_entidad) %>%
  summarise(y_h = weighted.mean(IRI, peso_aj), .groups="drop") %>%
  left_join(afijacion %>% select(departamento_entidad, W_h), by = "departamento_entidad") %>%
  summarise(IRI = sum(W_h * y_h)) %>% pull(IRI)

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("AJUSTE POR NO RESPUESTA\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("• Tasa respuesta:", round(nrow(muestra_resp)/nrow(muestra)*100, 1), "%\n")
cat("• IRI sin ajuste:", round(IRI_global, 4), "\n")
cat("• IRI ajustado:  ", round(IRI_aj, 4), "\n")
cat("• Diferencia:    ", round(abs(IRI_aj-IRI_global)/IRI_global*100, 2), "%\n")
cat("• Conclusión:     Impacto", ifelse(abs(IRI_aj-IRI_global)/IRI_global*100 < 2, "DESPRECIABLE", "MODERADO"), "\n\n")

# EXPORTAR
resultados <- list(
  muestra = muestra,
  estimacion_global = tibble(IRI = IRI_global, SE, CV),
  estimacion_dominios = est_dominios,
  dominios_constantes = dominios_constantes,
  ajuste_nr = ajuste
)

cat("═══════════════════════════════════════════════════════════════\n")
cat("✓ Resultados exportados en objeto 'resultados'\n")
cat("  - resultados$estimacion_dominios: Incluye CV corregido\n")
cat("  - resultados$dominios_constantes: Dominios con IRI homogéneo\n")
cat("═══════════════════════════════════════════════════════════════\n")

install.packages("kableExtra")
library("kableExtra")

## tabla 1
tabla_comp <- tibble(
  Diseño = c("MAS", "Estratificado"),
  Cobertura_Departamentos = c(
    n_distinct(muestra_mas$departamento_entidad),
    n_distinct(muestra$departamento_entidad)
  ),
  Dominios_n_ge_5 = c(dom_mas, dom_est),
  DEFF = c(1.0, 0.85)
)

kable(tabla_comp, caption = "Comparación de Diseños Muestrales",
      digits = 2) %>%
  kable_classic(full_width = F)
ggplot(tabla_comp, aes(x=Diseño, y=Dominios_n_ge_5, fill=Diseño)) +
  geom_col() +
  geom_text(aes(label=Dominios_n_ge_5), vjust=-0.3) +
  labs(title="Comparación de dominios estimables",
       y="Número de dominios estimables (n ≥ 5)",
       x=NULL) +
  theme_minimal()
## tabla 2
tabla_global <- tibble(
  Indicador = c("IRI Nacional", "Error Estándar", "Coeficiente de Variación (%)",
                "IC 95% (Límite Inferior)", "IC 95% (Límite Superior)"),
  Valor = c(
    round(IRI_global, 4),
    round(SE, 5),
    round(CV, 2),
    round(IRI_global - 1.96*SE, 4),
    round(IRI_global + 1.96*SE, 4)
  )
)

kable(tabla_global, caption="Estimación Nacional del IRI") %>%
  kable_classic(full_width=F)
ggplot(muestra, aes(x=IRI)) +
  geom_histogram(bins=30, fill="steelblue", color="white", alpha=0.8) +
  labs(title="Distribución del Índice de Riesgo (IRI) en la muestra",
       x="IRI", y="Frecuencia") +
  theme_minimal()

#gttabla 3
kable(top_dominios, caption="Top 15 Dominios con Mayor Riesgo (IRI)") %>%
  kable_classic(full_width=F)
ggplot(top_dominios, 
       aes(x=reorder(paste(Depto, Modal, sep=" - "), IRI), y=IRI)) +
  geom_col(fill="firebrick") +
  coord_flip() +
  labs(title="Top 15 dominios con mayor IRI",
       x="Dominio (Depto - Modalidad)",
       y="IRI promedio") +
  theme_minimal()


##tabla 4

ggplot(est_dominios, 
       aes(x=modalidad_de_contratacion, 
           y=departamento_entidad, fill=media_IRI)) +
  geom_tile(color="white") +
  scale_fill_gradient(low="white", high="red") +
  labs(title="Mapa de calor del IRI por Departamento × Modalidad",
       x="Modalidad de contratación", y="Departamento", fill="IRI") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))
#tabla 5

ggplot(est_dominios, aes(x=media_IRI, y=CV_final)) +
  geom_point(alpha=0.6, color="darkblue") +
  geom_smooth(method="lm", color="red", se=FALSE) +
  labs(title="Relación entre IRI y precisión (CV%)",
       x="IRI promedio del dominio",
       y="CV (%)") +
  theme_minimal()
#tabla 6


kable(table_precision, caption="Distribución de tipos de precisión por dominio") %>%
  kable_classic(full_width = F)

ggplot(est_dominios, aes(x=precision)) +
  geom_bar(fill="steelblue") +
  labs(title="Clasificación de precisión por dominios",
       x="Categoría de precisión",
       y="Número de dominios") +
  theme_minimal()

#tabla 7
ggplot(muestra, aes(x=departamento_entidad, y=IRI)) +
  geom_boxplot(fill="orange", alpha=0.7) +
  coord_flip() +
  labs(title="Distribución del IRI por Departamento",
       x="Departamento", y="IRI") +
  theme_minimal()
ggplot(muestra, aes(x=modalidad_de_contratacion, y=IRI)) +
  geom_boxplot(fill="skyblue") +
  coord_flip() +
  labs(title="Distribución del IRI por Modalidad",
       x="Modalidad", y="IRI") +
  theme_minimal()


library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(dplyr)

colombia <- ne_states(country = "Colombia", returnclass = "sf")

normalize <- function(x) {
  x %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    toupper() %>%
    trimws()
}

colombia$NAME_1_norm <- normalize(colombia$name)

df_mapa <- est_estratos %>%
  select(departamento_entidad, y_h) %>%
  mutate(depto_norm = normalize(departamento_entidad))

mapa_iri <- colombia %>%
  left_join(df_mapa, by = c("NAME_1_norm" = "depto_norm"))

ggplot(mapa_iri) +
  geom_sf(aes(fill = y_h), color = "gray20", size = 0.3) +
  scale_fill_gradient(
    low = "#FEE5D9",
    high = "#A50F15",
    name = "IRI promedio"
  ) +
  labs(
    title = "Índice de Riesgo IRI por Departamento",
    subtitle = "Estimación basada en muestreo estratificado",
    caption = "Fuente: SECOP 2020–2024 | Estimación propia"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  )