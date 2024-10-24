#let authorInfo = "Marcin Bajkowski, 193696"
#let dateInfo = datetime(year: 2024, month: 04, day: 28)
#let dateString = "28.04.2024r."

#set document(
  title: "Układy równań liniowych",
  author: authorInfo,
  date: dateInfo,
)

#set page(
  paper: "a4",
  margin: (x: 1.5cm, y: 1.5cm),
  footer: context [
    *#authorInfo, #dateString*
    #h(1fr)
    #counter(page).display(
      "1 / 1",
      both: true,
    )  
  ]
)

#align(center)[
  #stack(
    dir: ttb,
    text(size: 20pt)[*Układy równań liniowych*],
    v(10pt),
    text(size: 16pt)[Projekt 2 - Metody Numeryczne],
    v(10pt),
    text(size: 12pt)[#authorInfo, #dateString]
  )
]

= 1. Wstęp teoretyczny

Celem projektu jest implementacja i analiza dwóch metod iteracyjnych (Jacobiego i Gaussa-Seidla) oraz jednej metody bezpośredniej (faktoryzacja LU) rozwiązywania układów równań liniowych. Testy poszczególnych metod będą przeprowadzane na macierzach o różnych rozmiarach, które mogą powstać w wyniku dyskretyzacji równań różniczkowych i są powszczechnie stosowane w takich zagadnieniach jak: elektronika, elektrodynamika, mechanika (zastosowania lotnicze, biomechanika, motoryzacja), badanie wytrzymałości materiałów i konstrukcji, symulacje odkształceń, naprężeń, przemieszczeń i drgań, akustyka, fotonika, termodynamika, dynamika płynów i wiele innych.

Metody te zostały zaimplementowane w języku Python z wykorzystaniem własnoręcznie napisanej klasy *`Matrix`* reprezentującej macierz oraz bibliotek: *`time`* do pomiaru czasu wykonywania się poszczególnych metod, *`matplotlib`* do wizualizacji wyników, *`math`* do obliczeń matematycznych, *`os`* do operacji na plikach i folderach oraz *`copy`* do kopiowania obiektów. W dalszej części sprawozdania zostaną przedstawione wyniki przeprowadzonych testów wydajnościowych, które pozwolą na porównanie zaimplementowanych metod pod względem szybkości oraz dokładności rozwiązywania układów równań liniowych. Dodatkowo każda z metod została przetestowana na macierzach o różnych rozmiarach, aby sprawdzić jak zachowuje się w zależności od wielkości danych wejściowych. 

= 2. Wstęp praktyczny (metody)

Wszystkie wspomniane metody służą do rozwiązywania układów równań liniowych postaci $A x = b$, gdzie $A$ jest macierzą współczynników, $x$ wektorem niewiadomych, a $b$ wektorem wyrazów wolnych. W przypadku metody Jacobiego, rozwiązanie układu równań liniowych polega na iteracyjnym przybliżaniu wartości wektora $x$ na podstawie rekurencyjnego równania:
$
x #sub[i] #super[(k+1)] = 1/(a #sub[ii]) (b #sub[i] - sum_(j!=i) a #sub[ij] * x #sub[j] #super[(k)]), i = 1, 2, ..., n
$
Algorytm kończy się, gdy spełniony jest warunek zbieżności, czyli gdy norma błędu jest mniejsza od zadanego epsilon. Metoda ta nie jest jednak zawsze zbieżna, więc nie gwarantuje znalezienia rozwiązania. W takim przypadku należy zastosować inną metodę iteracyjną, np. Gaussa-Seidla. W tej metodzie wartość wektora $x$ jest obliczana na podstawie równania:
$
x #sub[i] #super[(k+1)] = 1/(a #sub[ii]) (b #sub[i] - sum_(j < i) a #sub[ij] * x #sub[j] #super[(k+1)] - sum_(j > i) a #sub[ij] * x #sub[j] #super[(k)]), i = 1, 2, ..., n
$
Metoda Gaussa-Seidla jest zbieżna dla macierzy symetrycznych i dodatnio określonych. W przeciwnym przypadku może nie zwrócić poprawnego wyniku. Warto zauważyć, że metoda ta jest mniej kosztowna obliczeniowo niż metoda Jacobiego, ponieważ w każdej iteracji korzysta z poprawionych wartości wektora $x$.

Ostatnią rozpatrzaną metodą jest faktoryzacja LU, która polega na zdekomponowaniu macierzy $A$ na iloczyn dwóch macierzy trójkątnych: dolnej $L$ i górnej $U$. Następnie rozwiązanie układu równań liniowych sprowadza się do rozwiązania dwóch układów równań liniowych z macierzami trójkątnymi. Metoda ta jest dokładna, ale wymaga więcej pamięci i czasu obliczeniowego niż metody iteracyjne. Do rozwiązania wykorzystuje się dwa wzory rekurencyjne:
$
L #sub[ii] * y #sub[i] = b #sub[i] - sum_(j=1) #super[(i-1)] L #sub[ij] * y #sub[j], i = 1, 2, ..., n
$
$
U #sub[ii] * x #sub[i] = y #sub[i] - sum_(j=i+1) #super[n] U #sub[ij] * x #sub[j], i = n, n-1, ..., 1
$

#pagebreak()

= 3. Dane wejściowe pierwszego zestawu danych

Zgodnie z treścią #underline[Zadania A] zaimplementowana została macierz $A$ o wymiarach $N #sym.times N$ dla indeksu $193696$.
Dane pozwalające wyliczyć konieczne wartości dla macierzy $A$ oraz wektora wyrazów wolnych $b$:
- $c$ - przedostatnia cyfra numeru indeksu ($c = 9$),
- $d$ - ostatnia cyfra numeru indeksu ($d = 6$),
- $e$ - czwarta cyfra numeru indeksu ($e = 6$),
- $f$ - trzecia cyfra numeru indeksu ($f = 3$),
- $a 1 = 5 + e$ - współczynniki macierzy $A$ ($a 1 = 11$),
- $a 2 = a 3 = $ -$1$ - współczynniki macierzy $A$,
Rozmiar macierzy definiowany jest wzorem $N = 9 c d$ co w przypadku indeksu $193696$ daje $N = 996$. Z kolei wektor wyrazów wolnych $b$ o długości $N$ dla $n$-tego elementu jest równy $sin(n #sym.dot (f + 1)) = sin(n)$.

$ A #sub[996 #sym.times 996] = mat(11, -1, -1, 0, 0, 0, 0, dots.h, 0, 0, 0, 0, 0, 0, 0;
          -1, 11, -1, -1, 0, 0, 0, dots.h, 0, 0, 0, 0, 0, 0, 0;
          -1, -1, 11, -1, -1, 0, 0, dots.h, 0, 0, 0, 0, 0, 0, 0;
          0, -1, -1, 11, -1, -1, 0, dots.h, 0, 0, 0, 0, 0, 0, 0;
          0, 0, -1, -1, 11, -1, -1, dots.h, 0, 0, 0, 0, 0, 0, 0;
          0, 0, 0, -1, -1, 11, -1, dots.h, 0, 0, 0, 0, 0, 0, 0;
          0, 0, 0, 0, -1, -1, 11, dots.h, 0, 0, 0, 0, 0, 0, 0;
          dots.v, dots.v, dots.v, dots.v, dots.v, dots.v, dots.v, dots.down, dots.v, dots.v, dots.v, dots.v, dots.v, dots.v, dots.v;
          0, 0, 0, 0, 0, 0, 0, dots.h, 11, -1, -1, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0, 0, dots.h, -1, 11, -1, -1, 0, 0, 0;
          0, 0, 0, 0, 0, 0, 0, dots.h, -1, -1, 11, -1, -1, 0, 0;
          0, 0, 0, 0, 0, 0, 0, dots.h, 0, -1, -1, 11, -1, -1, 0;
          0, 0, 0, 0, 0, 0, 0, dots.h, 0, 0, -1, -1, 11, -1, -1;
          0, 0, 0, 0, 0, 0, 0, dots.h, 0, 0, 0, -1, -1, 11, -1;
          0, 0, 0, 0, 0, 0, 0, dots.h, 0, 0, 0, 0, -1, -1, 11;
          )
  op(#h(10pt))
  b #sub[996]= mat(-0.7568;
          0.9894;
          -0.5366;
          -0.2879;
          0.9129;
          -0.9056;
          0.2709;
          dots.v;
          0.9997;
          -0.6365;
          -0.1677;
          0.8557;
          -0.951;
          0.3875;
          0.4444;
          )         
$

= 4. Rezultaty metod rozwiązywania układów równań liniowych dla macierzy zdefiniowanej w #underline[zadaniu A]

#underline[Zadania B] polegało na zaimplementowaniu metod iteracyjnych rozwiązywania układów równań liniowych: Jacobiego i Gaussa-Seidla (uzupełnione o #underline[Zadanie D]) oraz porównanie ich efektywności pod względem liczby iteracji potrzebnych do uzyskania rozwiązania o normie residuum mniejszej niż zadanego epsilon. Wartość epsilon została ustalona na $10#super[-9]$. Wyniki prezentują się następująco:

#align(center)[
  #figure(
    table(
      columns: 4,
      table.header(
        [*Metoda*], [*Liczba iteracji*], [*Czas obliczeń [s]*], [*Norma residuum*]
      ),
      text[*Jacobi*], text[$25$], text[$6.7150$], text[$8.278375316734148#sym.times 10#super[-10]$],
      text[*Gauss-Seidel*], text[$17$], text[$3.8724$], text[$3.8123459240370307#sym.times 10#super[-10]$],
      text[*Faktoryzacja LU*], text[$-$], text[$17.6411$], text[$2.34135805584451#sym.times 10#super[-15]$]
    ),
    caption: [Porównanie norm residuum dla metod iteracyjnych i bezpośredniej]
  )
]

Wyniki testów pokazują, że metoda Gaussa-Seidla jest szybsza ($3.8724s$) i dokładniejsza ($17$ iteracji) od metody Jacobiego ($6.7150s$ dla $25$ iteracji), co jest zgodne z oczekiwaniami. Róznicę w czasie obliczeń wynosi $2.8426s$, co jest znaczącą różnicą dla tak krótkich czasów. Z kolei metoda faktoryzacji LU okazała się najdokładniejsza ($2.34135805584451#sym.times 10#super[-15]$), ale jednocześnie najwolniejsza ($17.6411s$), co wynika z konieczności faktoryzacji macierzy $A$ na macierze $L$ i $U$.

Na podstawie tabeli można stwierdzić, że metoda Gaussa-Seidla jest najlepszym wyborem dla rozwiązywania układów równań liniowych, ponieważ jest zarówno szybsza, jak i dokładniejsza od metody Jacobiego. Metoda faktoryzacji LU jest dokładna, ale kosztowna obliczeniowo, dlatego nie jest zalecana dla dużych macierzy.

#pagebreak()

#figure(
  image("images/exercise_B.png"),
  caption: [Porównanie norm residuum dla metod iteracyjnych],
)

Wykres przedstawia porównanie norm residuum dla metod iteracyjnych: Jacobiego i Gaussa-Seidla. Jak widać, metoda Gaussa-Seidla osiąga mniejszą normę residuum w krótszym czasie niż metoda Jacobiego, co potwierdza wyniki z tabeli.

= 5. Dane wejściowe drugiego zestawu danych
Czy metody iteracyjne dla takich wartości elementów macierzy
A zbiegają się? Dla obu metod przedstaw na wykresie jak zmienia się norma residuum w
kolejnych iteracjach (oś y w skali logarytmicznej).

#underline[Zadanie C] polegało na stworzeniu układu równań z nowymi wartościami współczynników macierzy $A$ ($a 1 = 3$, $a 2 = a 3 = -1$) oraz obliczeniu wektora wyrazów wolnych $b$ zgodnie z treścią #underline[Zadania A].

$ A #sub[996 #sym.times 996] = mat(3, -1, -1, 0, 0, 0, 0, dots.h, 0, 0, 0, 0, 0, 0, 0;
          -1, 3, -1, -1, 0, 0, 0, dots.h, 0, 0, 0, 0, 0, 0, 0;
          -1, -1, 3, -1, -1, 0, 0, dots.h, 0, 0, 0, 0, 0, 0, 0;
          0, -1, -1, 3, -1, -1, 0, dots.h, 0, 0, 0, 0, 0, 0, 0;
          0, 0, -1, -1, 3, -1, -1, dots.h, 0, 0, 0, 0, 0, 0, 0;
          0, 0, 0, -1, -1, 3, -1, dots.h, 0, 0, 0, 0, 0, 0, 0;
          0, 0, 0, 0, -1, -1, 3, dots.h, 0, 0, 0, 0, 0, 0, 0;
          dots.v, dots.v, dots.v, dots.v, dots.v, dots.v, dots.v, dots.down, dots.v, dots.v, dots.v, dots.v, dots.v, dots.v, dots.v;
          0, 0, 0, 0, 0, 0, 0, dots.h, 3, -1, -1, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0, 0, dots.h, -1, 3, -1, -1, 0, 0, 0;
          0, 0, 0, 0, 0, 0, 0, dots.h, -1, -1, 3, -1, -1, 0, 0;
          0, 0, 0, 0, 0, 0, 0, dots.h, 0, -1, -1, 3, -1, -1, 0;
          0, 0, 0, 0, 0, 0, 0, dots.h, 0, 0, -1, -1, 3, -1, -1;
          0, 0, 0, 0, 0, 0, 0, dots.h, 0, 0, 0, -1, -1, 3, -1;
          0, 0, 0, 0, 0, 0, 0, dots.h, 0, 0, 0, 0, -1, -1, 3;
          )
  op(#h(10pt))
  b #sub[996]= mat(-0.7568;
          0.9894;
          -0.5366;
          -0.2879;
          0.9129;
          -0.9056;
          0.2709;
          dots.v;
          0.9997;
          -0.6365;
          -0.1677;
          0.8557;
          -0.951;
          0.3875;
          0.4444;
          )         
$

= 6. Rezultaty metod rozwiązywania układów równań liniowych dla macierzy zdefiniowanej w #underline[zadaniu C]

#underline[Zadanie C] polegało na zastosowaniu wcześniej zaimplementowanych metod iteracyjnych do rozwiązania układu równań liniowych z nowymi wartościami współczynników macierzy $A$ ($a 1 = 3$, $a 2 = a 3 = -1$) oraz obliczeniu wektora wyrazów wolnych $b$ zgodnie z treścią #underline[Zadania A]. Zadanie zostało uzupełnione dodatkowo o wyniki dla metody faktoryzacji LU (#underline[Zadanie D]). Poniżej przedstawiono porównanie ich efektywności pod względem liczby iteracji potrzebnych do uzyskania rozwiązania o normie residuum mniejszej niż zadanego epsilon oraz czasu obliczeń. Wartość epsilon pozostała niezmienna ($10#super[-9]$), a limit błędu został ustalony na $10#super[9]$. 

#pagebreak()

Wyniki po wprowadzonych zmianach prezentują się następująco:

#align(center)[
  #figure(
    table(
      columns: 5,
      table.header(
        [*Metoda*], [*Liczba iteracji*], [*Czas obliczeń [s]*], [*Norma residuum*], [*Zbieżność*]
      ),
      text[*Jacobi*], text[$60$], text[$15.9803$], text[$1301276180.7303445$], text[$times$],
      text[*Gauss-Seidel*], text[$24$], text[$5.8992$], text[$1029036099.3947304$], text[$times$],
      text[*Faktoryzacja LU*], text[$-$], text[$17.6019$], text[$1.9769795577221488#sym.times 10#super[-13]$], text[$-$]
    ),
    caption: [Porównanie norm residuum dla metod iteracyjnych i bezpośredniej]
  )
]

Wyniki testów pokazują, że metoda Gaussa-Seidla jest szybsza ($5.8992s$) i dokładniejsza ($24$ iteracje) od metody Jacobiego ($15.9803s$ dla $60$ iteracji), co jest zgodne z oczekiwaniami. Różnica w czasie obliczeń wynosi $9.0811s$, co jest ogromną różnicą dla tak krótkich czasów. Metoda faktoryzacji LU okazała się najdokładniejsza ($1.9769795577221488#sym.times 10#super[-13]$), ale jednocześnie najwolniejsza ($17.6019s$), co wynika z konieczności faktoryzacji macierzy $A$ na macierze $L$ i $U$. Co więcej, Jacobi okazał się być szybszy od faktoryzacji LU o tylko $1.3784s$, co jest niewielką różnicą w porównaniu do Gaussa-Seidla.

Warty uwagi jest fakt, iż Jacobi i Gauss-Seidel nie zbiegają się dla nowych wartości współczynników macierzy $A$ w porównaniu do poprzednich wartości zdefiniowanych w #underline[Zadaniu A].

#figure(
  image("images/exercise_C.png"),
  caption: [Porównanie norm residuum dla metod iteracyjnych],
)

Wykres przedstawia porównanie norm residuum dla metod iteracyjnych: Jacobiego i Gaussa-Seidla. Jak widać, metoda Gaussa-Seidla osiąga mniejszą normę residuum w krótszym czasie niż metoda Jacobiego, co potwierdza wyniki z tabeli.

= 7. Porównanie czasów obliczeń dla metod rozwiazywania układów równań liniowych

#underline[Zadanie E] polegało na utworzeniu wykresu prezentującego, jak zmienia się czas potrzebny na wyznaczenie rozwiązania dla trzech badanych metod (metoda Jacobiego, metoda Gaussa-Seidla oraz faktoryzacja LU) w zależności od liczby niewiadomych $N$. Przyjęto, że wartości $N$ miały przyjmować kolejne wartości z sekwencji ${100, 500, 1000, 2000, 3000}$, które odzwierciedlały rosnący rozmiar macierzy opisanej w #underline[Zadaniu A]. 

To zadanie pozwala na zrozumienie, jak skalowanie rozmiaru macierzy wpływa na czas wykonania każdej z badanych metod, co jest istotne w kontekście optymalizacji wydajności obliczeniowej. Co więcej pozwala na przyszły dobór odpowiedniej metody w zależności od rozmiaru danych wejściowych, problemu do rozwiązania, dostępnych zasobów obliczeniowych, wymagań dotyczących dokładności rozwiązania oraz czasu obliczeń.

#pagebreak()

Wyniki dla poszczególnych metod przedstawiono w tabeli poniżej:

#align(center)[
  #figure(
    table(
      columns: 4,
      table.header(
        [*Rozmiar macierzy*], [*Jacobi [s]*], [*Gauss-Seidel [s]*], [*Faktoryzacja LU [s]*]
      ),
      text[$bold(100)$], text[$0.0617$], text[$0.0392$], text[$0.0181$],
      text[$bold(500)$], text[$1.5391$], text[$0.9294$], text[$2.1670$],
      text[$bold(1000)$], text[$6.6282$], text[$3.8385$], text[$17.7077$],
      text[$bold(1500)$], text[$15.3376$], text[$9.2717$], text[$60.4937$],
      text[$bold(2000)$], text[$27.6398$], text[$16.6547$], text[$144.8885$],
      text[$bold(2500)$], text[$43.4565$], text[$26.2444$], text[$281.8669$],
      text[$bold(3000)$], text[$62.7243$], text[$37.5745$], text[$485.6950$],
    ),
    caption: [Porównanie czasów obliczeń dla metod rozwiązywania układów równań liniowych]
  )
]

#figure(
  image("images/exercise_E.png"),
  caption: [Porównanie czasów obliczeń dla metod rozwiązywania układów równań liniowych],
)

Wyniki testów pokazują, że metoda Gaussa-Seidla jest najszybsza dla wszystkich rozmiarów macierzy, co jest zgodne z wynikami z poprzednich testów. Metoda Jacobiego jest wolniejsza od metody Gaussa-Seidla, ale szybsza od faktoryzacji LU. W przypadku macierzy o rozmiarze $100 #sym.times 100$ czasy obliczeń wynoszą odpowiednio: $0.0617s$ dla metody Jacobiego, $0.0392s$ dla metody Gaussa-Seidla i $0.0181s$ dla metody faktoryzacji LU. Dla macierzy o rozmiarze $3000 #sym.times 3000$ czasy obliczeń wynoszą odpowiednio: $62.7243s$ dla metody Jacobiego, $37.5745s$ dla metody Gaussa-Seidla i $485.6950s$ dla metody faktoryzacji LU. Porównując czas potrzebny na wyznaczenie rozwiązania dla najmniejszej, jak i największej macierzy, można zauważyć, że czas obliczeń rośnie wraz z rozmiarem macierzy dla każdej z metod, co jest zgodne z oczekiwaniami. Róznice w czasie obliczeń dla najmniejszej macierzy wynosi $0.0211s$ na korzyść metody faktoryzacji LU w porównaniu z metodą Gaussa-Seidla, co jest niewielką różnicą. Natomiast dla największej macierzy różnica w czasie obliczeń wynosi $448.1205s$ na korzyść metody Gaussa-Seidla w porównaniu z metodą faktoryzacji LU, co jest znaczącą różnicą dla tak dużych czasów. Porównując metody Jacobiego i Gaussa-Seidla można zauważyć, że czas obliczeń rośnie liniowo wraz z rozmiarem macierzy i jest on większy średnio o $65%$ w stosunku do metody Jacobiego. Natomiast metoda faktoryzacji LU ma złożoność obliczeniową $O(n#super[3])$, co sprawia, że czas obliczeń rośnie kwadratowo wraz z rozmiarem macierzy.

= 8. Podsumowanie

Celem projektu było zaimplementowanie i analiza dwóch metod iteracyjnych (Jacobiego i Gaussa-Seidla) oraz jednej metody bezpośredniej (faktoryzacja LU) rozwiązywania układów równań liniowych. Testy poszczególnych metod zostały przeprowadzone na macierzach o różnych rozmiarach, aby porównać ich wydajność pod względem szybkości oraz dokładności rozwiązywania układów równań liniowych. 

Wyniki testów pozwoliły na stwierdzenie, że metoda Gaussa-Seidla jest najszybsza i najdokładniejsza dla rozwiązywania układów równań liniowych. Metoda Jacobiego jest wolniejsza od metody Gaussa-Seidla, ale szybsza od faktoryzacji LU. Metoda faktoryzacji LU jest najdokładniejsza, ale najwolniejsza, co wynika z konieczności faktoryzacji macierzy $A$ na macierze $L$ i $U$. 

W zależności od rozmiaru macierzy i wymagań dotyczących dokładności rozwiązania można wybrać odpowiednią metodę, aby zoptymalizować czas obliczeń. W przypadku metody Jacobiego i Gaussa-Seidla czas obliczeń rośnie liniowo wraz z rozmiarem macierzy, co jest związane z koniecznością wykonania większej liczby iteracji. Natomiast metoda faktoryzacji LU ma złożoność obliczeniową $O(n#super[3])$, co sprawia, że czas obliczeń rośnie kwadratowo wraz z rozmiarem macierzy. Porównując metody, można stwierdzić, że metoda Gaussa-Seidla jest najlepszym wyborem dla rozwiązywania układów równań liniowych, ponieważ jest zarówno szybsza, jak i dokładniejsza od metody Jacobiego. Metoda faktoryzacji LU jest dokładna, ale kosztowna obliczeniowo, dlatego nie jest zalecana dla dużych macierzy. W zależności od rozmiaru macierzy i wymagań dotyczących dokładności rozwiązania można wybrać odpowiednią metodę, aby zoptymalizować czas obliczeń.

Warto także zauważyć, że metody iteracyjne nie zawsze zbiegają się dla wszystkich wartości współczynników macierzy, co może prowadzić do braku rozwiązania. 

W przypadku macierzy o ogromnych rozmiarach, takich jak $3000 #sym.times 3000$, metoda Gaussa-Seidla jest najszybsza i w porównaniu z metodą faktoryzacji LU osiąga znacząco lepsze wyniki. Róznice te mogą być nawet kilkunastukrotne w zależności od rozmiaru macierzy. Dla testowanych danych największa różnica wyniosła $448.1205s$ na korzyść metody Gaussa-Seidla w porównaniu z metodą faktoryzacji LU (tym samym metoda faktoryzacji LU była wolniejsza o $1192,5877%$).

Warto zauważyć, że metody zostały napisane w języku Python z wykorzystaniem własnoręcznie napisanej klasy *`Matrix`* reprezentującej macierz. Wykorzystanie gotowych bibliotek takich jak *`numpy`* mogłoby przyspieszyć obliczenia, ale nie pozwoliłoby na pełne zrozumienie działania algorytmów. Innym sposobem na przyspieszenie obliczeń byłoby przepisanie kodu na przykład do języka C++, który jest znacznie szybszy od Pythona, ponieważ jest językiem kompilowanym, a nie interpretowanym. Dodatkowo każda funkcja w Pythonie jest obiektem, co wprowadza dodatkowe narzuty czasowe, które mogą wpłynąć na czas obliczeń.