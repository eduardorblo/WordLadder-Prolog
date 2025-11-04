% wordladder.pl
% Coloque 'lexico.txt' (cada palavra em uma linha) na mesma pasta.
% Execute em SWI-Prolog: ?- [wordladder].  then run main. or call main.

:- use_module(library(readutil)).   % read_line_to_string/2
:- use_module(library(apply)).      % maplist, foldl
:- use_module(library(assoc)).      % assoc for dist/conn maps
:- use_module(library(pairs)).      % pairs_keys_values/3

% ---------- Configurações ----------
word_len(3).      % comprimento das palavras (mudar conforme necessário)
max_iters(500).   % limite de iterações BFS

% ---------- Conversão palavra <-> número (8 bits por letra) ----------
% make_number(+AtomUpper, -Number)
make_number(Atom, Number) :-
    atom_codes(Atom, Codes),
    word_len(WL),
    length(Codes, WL),                 % garante comprimento correto
    make_number_codes(Codes, 0, 1, Number).

make_number_codes([], Num, _Mult, Num).
make_number_codes([C|Cs], Acc, Mult, Number) :-
    % espera A..Z (65..90)
    ( C >= 65, C =< 90 ->
        Byte is C - 65
    ;  % se não for A..Z, falha
       !, fail
    ),
    Acc2 is Acc + Byte * Mult,
    Mult2 is Mult * 256,
    make_number_codes(Cs, Acc2, Mult2, Number).

% make_word(+Number, -Atom)
make_word(Number, Atom) :-
    word_len(WL),
    make_word_chars(Number, WL, Codes),
    atom_codes(Atom, Codes).

make_word_chars(_Num, 0, []) :- !.
make_word_chars(Num, N, [C|Cs]) :-
    Byte is Num /\ 255,
    C is Byte + 65,
    Num2 is Num >> 8,
    N1 is N - 1,
    make_word_chars(Num2, N1, Cs).

% ---------- Pair LUT ----------
build_pair_lut(PairLut) :-
    word_len(WL),
    MaxIndex is WL - 1,
    findall(Val,
            ( between(0, MaxIndex, I),
              between(1, 255, J),      % 0 não altera, por isso ignora
              Shift is J << (I * 8),
              Val = Shift
            ),
            L),
    sort(L, PairLut).

% ---------- Carregar dicionário (mantendo originais com acentos) ----------
% load_words(+FileBase, -WordsNums, -OrigAtoms)
% WordsNums: lista de números usados internamente
% OrigAtoms: lista de átomos originais (com acentos), alinhados com WordsNums
load_words(FileBase, WordsNums, OrigAtoms) :-
    atom_concat(FileBase, '.txt', FileName),
    catch(open(FileName, read, In, [type(text)]), E,
          ( format('ERRO abrindo ~w: ~w~n', [FileName,E]), fail )),
    set_stream(In, encoding(utf8)),
    read_all_lines(In, Lines),
    close(In),

    word_len(WL),
    % filtra linhas com comprimento correto (em número de caracteres)
    include(length_codes_eq(WL), Lines, LenLines),

    % trim e uppercase (mantendo acentos nas originais)
    maplist(string_trim, LenLines, Trimmed),
    maplist(string_upper, Trimmed, UpLines),

    % OrigAtomsCandidate: átomos originais (MAIÚSCULOS, com acentos)
    maplist(atom_string, OrigAtomsCandidate, UpLines),

    % NormLinesCandidate: converter cada UpLine para versão ASCII por tabela
    maplist(line_to_norm_maybe, UpLines, NormMaybe),

    % NormMaybe is list of either norm(String) or invalid
    % Build parallel lists of only valid pairs: NormAtom and OrigAtom
    pair_filter_make_numbers(OrigAtomsCandidate, NormMaybe, WL, PairsNumAtom),

    % ordenar por número (chave) e extrair listas alinhadas
    keysort(PairsNumAtom, SortedPairs),        % SortedPairs: Num-Orig pairs
    pairs_keys_values(SortedPairs, Nums0, Atoms0),

    % remover possíveis duplicatas numéricas (mantendo primeira ocorrência)
    remove_duplicate_numbers(Nums0, Atoms0, WordsNums, OrigAtoms).

% Lê todas as linhas de um stream para uma lista de strings
read_all_lines(Stream, Lines) :-
    read_line_to_string(Stream, S),
    ( S == end_of_file -> Lines = []
    ; Lines = [S|Rest], read_all_lines(Stream, Rest)
    ).

length_codes_eq(WL, S) :-
    string_codes(S, Codes),
    length(Codes, L),
    L =:= WL.

string_to_atom_trim(S, Atom) :-
    string_codes(S, Codes),
    % remove possíveis espaços laterais
    string_trim(S, ST),
    atom_string(Atom, ST).

string_trim(Str, Trim) :-
    string_codes(Str, Cs),
    trim_codes(Cs, CsTrim),
    string_codes(Trim, CsTrim).

trim_codes(Cs, Out) :-
    trim_leading(Cs, C1),
    trim_trailing(C1, Out).

trim_leading([C|Cs], Out) :- C =< 32, !, trim_leading(Cs, Out).
trim_leading(Cs, Cs).

trim_trailing(Cs, Out) :- reverse(Cs, Rc), trim_leading(Rc, R1), reverse(R1, Out).

atom_is_alpha_len(WL, Atom) :-
    atom_codes(Atom, Codes),
    length(Codes, WL),
    forall(member(C, Codes), (C >= 65, C =< 90)).

% ---------- Transliteration table (Latin-1 & common) ----------
% map_accent_code(+CodePoint, -BaseCode) where BaseCode is ASCII uppercase A..Z (65..90)
% This table maps many Latin-1 accented uppercase/lowercase letters to base A..Z.
map_accent_code(193, 65).  % Á -> A
map_accent_code(194, 65).  % Â -> A
map_accent_code(195, 65).  % Ã -> A
map_accent_code(199, 67).  % Ç -> C
map_accent_code(201, 69).  % É -> E
map_accent_code(202, 69).  % Ê -> E
map_accent_code(205, 73).  % Í -> I
map_accent_code(206, 73).  % Î -> I
map_accent_code(211, 79).  % Ó -> O
map_accent_code(212, 79).  % Ô -> O
map_accent_code(213, 79).  % Õ -> O
map_accent_code(218, 85).  % Ú -> U
map_accent_code(219, 85).  % Û -> U

% lowercase Latin-1

map_accent_code(225, 65).  % á -> A
map_accent_code(226, 65).  % â -> A
map_accent_code(227, 65).  % ã -> A
map_accent_code(231, 67).  % ç -> C
map_accent_code(233, 69).  % é -> E
map_accent_code(234, 69).  % ê -> E
map_accent_code(237, 73).  % í -> I
map_accent_code(238, 73).  % î -> I
map_accent_code(243, 79).  % ó -> O
map_accent_code(244, 79).  % ô -> O
map_accent_code(245, 79).  % õ -> O
map_accent_code(250, 85).  % ú -> U
map_accent_code(251, 85).  % û -> U


% ---------- Converte uma linha em versão "norm" ASCII (ou falha) ----------
% line_to_norm_maybe(+LineString, -NormStringOrInvalid)
% Se a linha puder ser totalmente transliterada para A..Z, retorna norm(OutString)
% Caso contrário retorna invalid
line_to_norm_maybe(Line, NormStr) :-
    string_codes(Line, Codes),
    try_map_codes(Codes, MappedCodes),   % falha se algum char não mapeável
    string_codes(NormStr, MappedCodes), !.
line_to_norm_maybe(_Line, invalid).

% try_map_codes(+Codes, -MappedCodes) succeeds only if all codes can be mapped to A..Z codes.
try_map_codes([], []).
try_map_codes([C|Cs], [M|Ms]) :-
    map_code_to_ascii(C, M), !,
    try_map_codes(Cs, Ms).

% map_code_to_ascii(+Code, -AsciiUpperCode)
% - If Code already A..Z, keep it.
% - If Code a..z, convert to upper.
% - Else if appears in map_accent_code/2, use mapped base.
map_code_to_ascii(Code, Code) :-
    Code >= 65, Code =< 90.   % 'A'..'Z'
map_code_to_ascii(Code, Up) :-
    Code >= 97, Code =< 122,  % 'a'..'z' -> uppercase
    Up is Code - 32.
map_code_to_ascii(Code, Base) :-
    map_accent_code(Code, Base).

% ---------- Helpers para filtragem e conversão ----------
% pair_filter_make_numbers(+OrigAtomsList, +NormMaybeList, +WL, -Pairs)
% Produz lista de pares Num-OrigAtom para NormMaybe entries válidos
pair_filter_make_numbers([], [], _WL, []).
pair_filter_make_numbers([Orig|Origs], [Norm|Norms], WL, Pairs) :-
    ( Norm \= invalid,
      string_codes(Norm, NormCodes),
      length(NormCodes, WL),
      atom_string(NormAtom, Norm),
      atom_is_alpha_len(WL, NormAtom),
      catch(make_number(NormAtom, Num), _, fail) ->
        Pairs = [Num-Orig | Rest]
    ; Pairs = Rest
    ),
    pair_filter_make_numbers(Origs, Norms, WL, Rest).

% remove_duplicate_numbers(+Nums, +Atoms, -UniqueNums, -UniqueAtoms)
remove_duplicate_numbers(Nums, Atoms, UniqueNums, UniqueAtoms) :-
    remove_duplicate_numbers_(Nums, Atoms, [], RevNums, [], RevAtoms),
    reverse(RevNums, UniqueNums),
    reverse(RevAtoms, UniqueAtoms).

remove_duplicate_numbers_([], [], _Seen, [], _SeenAtoms, []).
remove_duplicate_numbers_([N|Ns], [A|As], Seen, OutNums, SeenAtoms, OutAtoms) :-
    ( memberchk(N, Seen) ->
        remove_duplicate_numbers_(Ns, As, Seen, OutNums, SeenAtoms, OutAtoms)
    ; remove_duplicate_numbers_(Ns, As, [N|Seen], OutNums2, [A|SeenAtoms], OutAtoms2),
      OutNums = [N|OutNums2],
      OutAtoms = [A|OutAtoms2]
    ).

% ---------- Encontrar pares ----------
find_all_pairs(WordsNums, PairLut, Pairs) :-
    findall(pair(W1Atom, W2Atom),
            ( nth0(I, WordsNums, W1),
              nth0(J, WordsNums, W2),
              J < I,
              Xor is W1 xor W2,
              memberchk(Xor, PairLut),
              make_word(W1, W1Atom),
              make_word(W2, W2Atom)
            ),
            Pairs).

% ---------- Escrever DOT (usa OrigAtoms para rótulos se disponível) ----------
% write_dot(+WordsNums, +Pairs, +OutFile, +OrigAtoms)
write_dot(WordsNums, Pairs, OutFile, OrigAtoms) :-
    catch(open(OutFile, write, Out, [encoding(utf8)]), E,
          ( format('ERRO abrindo ~w para escrita: ~w~n', [OutFile,E]), fail )),
    write(Out, 'graph words {\n'),
    % escrever nós: usar OrigAtoms alinhado quando possível
    forall(nth0(I, WordsNums, W),
           (
               ( nth0(I, OrigAtoms, Orig) -> fix_keyword(Orig, Safe) ; make_word(W, Atom), fix_keyword(Atom, Safe) ),
               format(Out, '  "~w";~n', [Safe])
           )),
    % escrever arestas
    forall(member(pair(Aatom,Batom), Pairs),
           (
               ( catch(make_number(Aatom, NA), _, fail), nth0(IndexA, WordsNums, NA) ->
                   nth0(IndexA, OrigAtoms, OrigA), ( OrigA = _ -> fix_keyword(OrigA, SA) ; make_word(NA, MA), fix_keyword(MA, SA) )
               ; ( make_word_atom_fallback(Aatom, SA) )
               ),
               ( catch(make_number(Batom, NB), _, fail), nth0(IndexB, WordsNums, NB) ->
                   nth0(IndexB, OrigAtoms, OrigB), ( OrigB = _ -> fix_keyword(OrigB, SB) ; make_word(NB, MB), fix_keyword(MB, SB) )
               ; ( make_word_atom_fallback(Batom, SB) )
               ),
               format(Out, '  "~w" -- "~w";~n', [SA, SB])
           )),
    write(Out, '}\n'),
    close(Out).

% helper fallback for atoms that are not convertible (shouldn't usually happen)
make_word_atom_fallback(AtomIn, Safe) :-
    ( catch(make_number(AtomIn, N), _, fail) ->
        make_word(N, A), fix_keyword(A, Safe)
    ; fix_keyword(AtomIn, Safe)
    ).

keywords(['NODE','EDGE','GRAPH','DIGRAPH','SUBGRAPH','STRICT']).
fix_keyword(Atom, Fixed) :-
    ( keywords(Ks), memberchk(Atom, Ks) -> atom_concat('_', Atom, Fixed) ; Fixed = Atom ).

% ---------- BFS por camadas (iterativo) ----------
% compute_distances(+AllWords, +PairLut, +Target, +MaxIters, -DistAssoc, -ConnAssoc)
compute_distances(AllWords, PairLut, Target, MaxIters, FinalDistAssoc, FinalConnAssoc) :-
    empty_assoc(Empty),
    % inicializa dist com -1 para todas as palavras
    foldl(init_word_assoc, AllWords, Empty, Dist0),
    put_assoc(Target, Dist0, 0, Dist1),
    empty_assoc(EmptyC),
    bfs_layers([Target], 0, MaxIters, AllWords, PairLut, Dist1, EmptyC, FinalDistAssoc, FinalConnAssoc).

init_word_assoc(W, AssocIn, AssocOut) :-
    put_assoc(W, AssocIn, -1, AssocOut).

% bfs_layers(Frontier, Iter, MaxIters, AllWords, PairLut, DistIn, ConnIn, DistOut, ConnOut)
bfs_layers(_, Iter, MaxIters, _All, _PairLut, Dist, Conn, Dist, Conn) :-
    Iter > MaxIters, !.
bfs_layers([], _Iter, _Max, _All, _PairLut, Dist, Conn, Dist, Conn) :- !.
bfs_layers(Frontier, Iter, MaxIters, AllWords, PairLut, DistIn, ConnIn, DistOut, ConnOut) :-
    NextIter is Iter + 1,
    findall(NewNodes-Parent,
            ( member(Parent, Frontier),
              get_assoc(Parent, DistIn, DistParent), DistParent =:= Iter,
              % para cada Parent, procure vizinhos não visitados
              findall(Neighbor,
                      ( member(Neighbor, AllWords),
                        get_assoc(Neighbor, DistIn, Dn),
                        Dn == -1,
                        Xor is Parent xor Neighbor,
                        memberchk(Xor, PairLut)
                      ),
                      Neighbors),
              member(NeighborNode, Neighbors),
              NewNodes = NeighborNode
            ),
            PairsFlat),
    % PairsFlat é uma lista de New-Parent entries; precisamos processar unicamente
    ( PairsFlat == [] ->
        DistOut = DistIn, ConnOut = ConnIn
    ; process_pairs(PairsFlat, DistIn, ConnIn, Iter, DistNext, ConnNext),
      % prepare next frontier as unique nodes discovered at this layer
      assoc_to_list(DistNext, L),  % get list to filter
      findall(W, (member(W-D, L), D =:= NextIter), NextFrontier),
      bfs_layers(NextFrontier, NextIter, MaxIters, AllWords, PairLut, DistNext, ConnNext, DistOut, ConnOut)
    ).

% process_pairs(ListOfNew-ParentEntries, DistIn, ConnIn, Iter, DistOut, ConnOut)
process_pairs([], Dist, Conn, _Iter, Dist, Conn).
process_pairs([Node-Parent|Rest], DistIn, ConnIn, Iter, DistOut, ConnOut) :-
    ( get_assoc(Node, DistIn, Dn), Dn == -1 ->
        NewDist is Iter + 1,
        put_assoc(Node, DistIn, NewDist, DistTemp),
        put_assoc(Node, ConnIn, Parent, ConnTemp)
    ; DistTemp = DistIn, ConnTemp = ConnIn
    ),
    process_pairs(Rest, DistTemp, ConnTemp, Iter, DistOut, ConnOut).

% ---------- Impressão do caminho ----------
print_path(FromNum, ToNum, ConnAssoc, DistAssoc) :-
    ( get_assoc(FromNum, DistAssoc, Dist), Dist =:= -1 ->
        writeln('Can not connect!')
    ; format('~w steps~n', [Dist]),
      print_path_chain(FromNum, ToNum, ConnAssoc)
    ).

print_path_chain(Num, ToNum, ConnAssoc) :-
    make_word(Num, Atom),
    writeln(Atom),
    ( Num =:= ToNum -> true
    ; ( get_assoc(Num, ConnAssoc, Parent) ->
            print_path_chain(Parent, ToNum, ConnAssoc)
      ; writeln('Stuck (no parent)!')
      )
    ).

% ---------- Fluxo principal ----------
main :-
    word_len(WL),
    max_iters(MI),
    format('WORD_LEN = ~w~n', [WL]),
    format('MAX_ITERS = ~w~n', [MI]),
    build_pair_lut(PairLut),
    length(PairLut, L), format('Creating Diff Lookup Table (size ~w)...~n', [L]),

    format('Loading Dictionary...~n'),
    load_words('lexico', WordsNums, OrigAtoms),
    length(WordsNums, CountWords),
    format('Loaded ~w words.~n', [CountWords]),

    format('Finding All Connections...~n'),
    find_all_pairs(WordsNums, PairLut, Pairs),
    length(Pairs, PairCount),
    format('Found ~w connections.~n', [PairCount]),
    write_dot(WordsNums, Pairs, 'graph.dot', OrigAtoms),
    format('Wrote graph.dot~n~n'),

    interactive_loop(WordsNums, PairLut).

% ---------- Loop interativo ----------
interactive_loop(WordsNums, PairLut) :-
    repeat,
      write('From Word: '), flush_output, read_line_to_string(user_input, S1),
      ( S1 == end_of_file -> !, fail ; true ),
      normalize_input(S1, FromAtom),
      write('  To Word: '), flush_output, read_line_to_string(user_input, S2),
      ( S2 == end_of_file -> !, fail ; true ),
      normalize_input(S2, ToAtom),
      ( FromAtom == '' -> writeln('Empty from word, exiting.'), !, fail ; true ),
      ( catch(make_number(FromAtom, FromNum), _, fail) ->
          true
      ; format('Invalid From Word (wrong length or chars).~n'), fail ),
      ( catch(make_number(ToAtom, ToNum), _, fail) ->
          true
      ; format('Invalid To Word (wrong length or chars).~n'), fail ),

      ( \+ memberchk(FromNum, WordsNums) ->
          format('No connections to ~w~n', [FromAtom]), fail
      ; true ),
      ( \+ memberchk(ToNum, WordsNums) ->
          format('No connections to ~w~n', [ToAtom]), fail
      ; true ),

      max_iters(MaxIters),
      compute_distances(WordsNums, PairLut, ToNum, MaxIters, DistAssoc, ConnAssoc),
      ( get_assoc(FromNum, DistAssoc, DFrom), DFrom =:= -1 ->
            writeln('Can not connect!')
      ; print_path(FromNum, ToNum, ConnAssoc, DistAssoc)
      ),
    fail.

normalize_input(S, Atom) :-
    string_trim(S, T),
    ( T == "" -> Atom = '' ; string_upper(T, TU), atom_string(Atom, TU) ).
