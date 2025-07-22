# Análisis Bayesiano para la Detección de Firmas Espectrales de Fusarium en Plantas de Banano


El marchitamiento por Fusarium causado por Fusarium oxysporum f.sp. cubense representa una de las amenazas más devastadoras para la producción mundial de banano, con pérdidas económicas que pueden ser muy perjudiciales para las plantaciones afectadas. La detección temprana durante el período asintomático, que puede extenderse hasta 20 días post-inoculación, es crucial para implementar medidas de control efectivas antes de que los síntomas visibles se manifiesten. 


Los métodos tradicionales de detección, incluyendo inspección visual y técnicas moleculares como PCR y ELISA, presentan limitaciones significativas en términos de costo, tiempo y aplicabilidad en campo. La espectroscopía de reflectancia VIS/NIR ha emergido como una alternativa prometedora, ofreciendo detección no destructiva y en tiempo real de cambios fisiológicos asociados con infecciones patógenas (Díaz, 2024).
Trabajos previos han demostrado el potencial de la espectroscopía para la detección de enfermedades en plantas, pero típicamente se basan en enfoques frecuentistas que no incorporan conocimiento previo sobre la progresión temporal de las enfermedades. 


Este estudio presenta un enfoque bayesiano que busca mejorar la detección temprana de Fusarium en plantas de banano integrando información de la literatura sobre las etapas de desarrollo de Fusarium en modelos de selección de variables, logrando una caracterización espectral evolutiva de la enfermedad, identificando patrones temporales distintivos y ventanas críticas de detección durante el período asintomático. Con este fin, se ajusta un modelo de regresión logística Bayesiana con selección de variables para cada uno de los días post inoculación. La informaicón previa de la literatura se agrega al modelo en la configuración de los hiperparámetros de las distribuciones apriori. Posteriormente, se evalúa la calidad del enfoque por medio de la integración de los resultados con técnicas de aprendizaje automático avanzadas, entrenando un modelo XGBoost con los datos de los días 1 a 15 y otro modelo XGBoost con estos mismos datos, agregando carcaterísticas informadas por el análisis Bayesiano. 


Se utilizan los datos recolectados en (Díaz, 2024) donde se utilizaron 240 plantas de banano Gros Michel cultivadas bajo condiciones controladas en el invernadero de AUGURA, Carepa, Antioquia. Las plantas fueron sometidas a ocho tratamientos diferentes: control, estrés hídrico, FOCR1, Ralstonia solanacearum Raza 2, y sus interacciones. Los datos espectrales se obtuvieron mediante un espectrofotómetro portátil ASD FieldSpec 4 Hi-Res NG en el rango de 350-2500 nm durante 15 días post-inoculación.


El trabajo original* alcanza su mejor performance en el día 3, alcanzando un accuracy de 86%. El modelo Bayesiano alcanza mejores resultados en los días 3 y 6, alcanzando 91% de accurary en el día 6. Además, el modelo XGBoost con características informadas por el análisis Bayesiano logró un notable F1-score de 99,02%, mejorando el resultado del modelo sin caracteristicas.



*Trabajo original: Díaz Herrera, Cristian Camilo. Detección temprana de estrés biótico y abiótico usando modelos de clasificación de datos de espectroscopía de reflectancia VIS/NIR: Aplicación en plantas de Banano Gros Michel. Universidad Nacional de Colombia, Facultad de Ciencias, Departamento de Estadística. 2024.
