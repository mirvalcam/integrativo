HNSC

Set de datos de Head and Neck Squamous Cell Carcinoma (HNSC) con un total de 443 observaciones en múltiples ómicas.

Ómicas disponibles:

- Variables clínicas
- RNA-Seq
- miRNA
- CNV

- Variables clínicas
- RNA-Seq
- TF
- Mutaciones
- Archivo de asociones

Cada grupo dispone de un conjunto de ómicas con las que trabajar tal como se ha comentado en la primera clase. Además, se incluye información
sobre los genes y miRNA diferencialmente expresados.

Este dataset, en vez de estar dispuesto en un único fichero, las ómicas se encuentran por separado. Al cargar cada archivo
.RData en el entorno de R aparecerán dos data.frames denominados X e Y.

El data.frame X siempre será la ómica en cuestión, y el data.frame Y se trata de una matriz de supervivencia donde cada muestra
contiene dos variables, un tiempo que indica el tiempo que el sujeto ha estado en el estudio hasta que ha ocurrido el evento de 
interés o bien desde el inicio hasta el fin del estudio o pérdida de seguimiento. La segunda variable, indica si ese paciente llegó
a fallecer (valor TRUE) o no durante el estudio. La matrix Y de todos los RData es la misma, con una es suficiente.

La única excepción radica en el archivo de Factores de Transcripción donde se carga una única matrix denominada X_TF.

Dentro de las variables clínicas, la variable que necesita más explicación es la variable numérica "tobacco_smoking_history" 
que sigue la siguiente estructura:

1 = Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)
2 = Current smoker (includes daily smokers and non-daily smokers or occasional smokers)
3 = Current reformed smoker for > 15 years (greater than 15 years)
4 = Current reformed smoker for ≤15 years (less than or equal to 15 years)
5 = Current reformed smoker, duration not specified

Posibles objetivos del estudio multi-ómico:
-	Obtener una firma multi-ómica que diferencie entre los diferentes tipos de pacientes.
-	Realizar un análisis de supervivencia multi-ómico que explique la mortalidad de los pacientes.
-	Entender la regulación de la expresión génica mediada por CNVs y miRNAs en general o en distintos tipos de pacientes.
-	Entender la regulación de la expresión génica mediada por TF y mutaciones usando, si se puede, las asociaciones biológicas inscritas en el archivo de asociaciones.
