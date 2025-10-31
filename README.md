Word Ladder Solver (Prolog)

Este programa em Prolog, inspirado no código do CodeParade (https://github.com/HackerPoet/WordChainSolver), resolve o problema da Word Ladder usando um dicionário de palavras fornecido. O programa permite encontrar todas as sequências de palavras em que cada palavra difere da anterior por apenas uma letra.

Funcionalidades

- Permite especificar o tamanho da palavra (WORD_LEN) e número máximo de iterações (MAX_ITERS) para buscas.
- Gera um arquivo .dot opcional representando o grafo de conexões. (Use Gephi).
- Permite entrada interativa para descobrir caminhos entre palavras.

Requisitos

- SWI-Prolog (ou outro interpretador compatível com ISO Prolog)
  - Site oficial: https://www.swi-prolog.org
- Arquivo de dicionário de palavras (lexico.txt) no mesmo diretório
  - Cada palavra deve estar em uma linha.
  - Preferencialmente todas em maiúsculas ou o programa converte internamente.

Instruções de Execução

1. Abra o terminal no diretório do projeto.
2. Inicie o SWI-Prolog:

swipl

3. Carregue o programa:

?- [wordladder].

4. Execute a função principal:

?- main.

5. Siga as instruções interativas:

- Informe a palavra de origem (From Word)
- Informe a palavra de destino (To Word)
- O programa exibirá o caminho de palavras, onde cada palavra difere da anterior por uma letra.

Configurações

- WORD_LEN: Defina o tamanho das palavras que deseja analisar.  
- MAX_ITERS: Define a profundidade máxima da busca.  
- Arquivo de saída .dot: Pode ser usado no Gephi para visualizar o grafo.

Observações

- Para palavras longas ou grandes dicionários, o programa pode ser lento devido à complexidade de busca.  
- Recomenda-se iniciar com palavras de 3 ou 4 letras para testes.  

