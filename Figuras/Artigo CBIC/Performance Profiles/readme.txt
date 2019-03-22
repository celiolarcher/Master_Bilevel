
Ola,

Segue anexo uma implementação dos PPs (.jar) e um conjunto de arquivos de exemplo (.dat).
Compactei tudo em um arquivo ZIP.

Para utilizar a implementação, extraia os arquivos.


Para executar o exemplo, basta digitar o comando que segue no prompt (lembrando de estar na pasta dos arquivos extraidos).

java -jar Performance.Profiles.jar problems.dat algorithms.dat 1 example_pp.gnu -1


Para utilizar seus dados, o comando segue a regra seguinte:

java -jar Performance.Profiles.jar <arquivo.que.lista.os.problemas.dat> <arquivo.que.lista.os.algoritmos.dat> 1 <arquivo.de.saida.gnu> -1

onde
<arquivo.que.lista.os.problemas.dat> é o arquivo que contem os nomes dos problemas que serão incluídos na comparação.
<arquivo.que.lista.os.algoritmos.dat> é o arquivo que contem os nomes dos algoritmos que estão sendo comparados.
                                                       É importante perceber que para cada nome nesse arquivo deve existir um arquivo <algoritmo.dat> associado,
                                                       com os valores obtidos (seguindo a métrica em questão) para os problemas listados em <arquivo.que.lista.os.problemas.dat>.
                                                       Assim, tendo um algoritmo a0 e um problema p3, então estes devem estar nos respectivos arquivos (de algoritmos e problemas) e
                                                       deve-se ter um arquivo a0.dat com uma linha contendo o valor do desempenho obtido por esse algoritmo para o problema p3 na forma:
                                                       p3 valor.desempenho
                                                       por exemplo,
                                                       p3 234.567
                                                       Vide os arquivos do exemplo para facilitar o entendimento.
                                                       Note que cada linha deve conter o desempenho para cada problema e o separador do nome do problema e do valor do desempenho pode ser
                                                       espaço ou tabulação.
                                                       obs.: o nome do problema e do algoritmo não pode conter espaço.
<arquivo.de.saida.gnu> é o arquivo gnuplot que será gerado. Executando ele no gnuplot, 2 imagens (eps) são geradas: (i) com os performance profiles e (ii) com o gráfico de barras correspondente.

Note que essa implementação não gera os PPs, mas gera um arquivo gnuplot e os dados necessários para isso.
Também, apesar do exemplo já conter os dados normalizados, isso não é necessário (mas note que os valores de desempenho devem ser diferentes de zero e preferencialmente positivos).

Os números 1 e -1 da linha de execução são flags.
O primeiro número (1) define o tipo de PPs que se deseja gerar.
O valor 1 é o caso padrão (aconselho não mudar isso).
O segundo número (-1) define o valor máximo para o eixo x (ou seja, tau).
O valor -1 indica que o tau_{max} = máximo valor nos dados (ou seja, o maior valor possível).
Normalmente, eu deixo o -1 e mudo a "aparência" dos PPs no arquivo gnuplot, quando necessário (para melhorar algum aspecto visual).

Precisando de ajuda é só dizer.

obs.: aindei fazendo umas pequenas alterações na implementação para gerar uns dados extras (alguns deles, que acho interessantes, são impressos ao executar).
        retirei algumas coisas (p/ simplificar), mas acho que não gerei nenhum problema na implementação por conta disso.

Abraço,
Heder 
