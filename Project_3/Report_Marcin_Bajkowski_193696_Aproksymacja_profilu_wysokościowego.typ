#let authorInfo = "Marcin Bajkowski, 193696"
#let dateInfo = datetime(year: 2024, month: 05, day: 26)
#let dateString = "26.05.2024r."

#set document(
  title: "Aproksymacja profilu wysokościowego",
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
    text(size: 20pt)[*Aproksymacja profilu wysokościowego*],
    v(10pt),
    text(size: 16pt)[Projekt 3 - Metody Numeryczne],
    v(10pt),
    text(size: 12pt)[#authorInfo, #dateString]
  )
]

= 1. Wstęp teoretyczny

Celem projektu jest zaimplementowanie oraz porównanie algorytmów interpolacji funkcji na przykładzie aproksymacji profilu wysokościowego terenu. W ramach projektu zaimplementowano *Metodę interpolacji wielomianowej Lagrange’a*, *Metodę interpolacji wielomianowej zawierającej wezły Chebyshewa* oraz *Metodę interpolacji wykorzystującą funkcje sklejane trzeciego stopnia*. Testy poszczególnych metod zostały przeprowadzone na trzech różnych zestawach danych, aby zbadać zachowanie algorytmów w różnych warunkach. Ponadto, przeanalizowane zostaną różnice wynikające z różnych parametrów interpolacji, takich jak liczba punktów czy ich rozmieszczenie.

Metody te zostały zaimplementowane w języku Python z wykorzystaniem własnoręcznie napisanej klasy *`Matrix`* reprezentującej macierz oraz bibliotek: *`time`* do pomiaru czasu wykonywania się poszczególnych metod, *`matplotlib`* do wizualizacji wyników, *`math`* do obliczeń matematycznych, *`os`* do operacji na plikach i folderach oraz *`copy`* do kopiowania obiektów. W dalszej części sprawozdania zostaną przedstawione wyniki przeprowadzonych testów wydajnościowych, które pozwolą na porównanie zaimplementowanych metod pod względem szybkości oraz dokładności rozwiązywania układów równań liniowych. Dodatkowo każda z metod została przetestowana na macierzach o różnych rozmiarach, aby sprawdzić jak zachowuje się w zależności od wielkości danych wejściowych.

= 2. Wstęp praktyczny (metody)
== 2.1. Metoda interpolacji wielomianowej Lagrange’a (zawierająca wezły równoodległe oraz wezły Czebyszewa)

Wszystkie wspomniane metody służą do interpolacji funkcji na podstawie znanych wartości w punktach. Metoda interpolacji wielomianowej Lagrange’a polega na znalezieniu wielomianu stopnia $n$ przechodzącego przez $n + 1$ punktów $(x #sub[0], y #sub[0]), (x #sub[1], y #sub[1]), ..., (x #sub[n+1], y #sub[n+1])$. Bazą tej metody jest wzór *{1}*, a wielomian interpolacyjny Lagrange'a można zapisać jako sumę iloczynów wartości funkcji w punktach $y #sub[i]$ oraz funkcji bazowej $phi.alt_i (x)$ *{2}*. Wzory:
$ phi.alt_i (x) = limits(product)_(j=1,j!=i)^(n+1) frac((x - x_j), (x_i - x_j))
op(#h(50pt))
F(x) = sum_(i=1)^(n+1) y_i phi.alt_i (x) $

Metoda Lagrange’a ma również wady, takie jak podatność na efekt Rungego, który polega na oscylacjach wielomianu interpolacyjnego w okolicach krańców przedziału. W celu zminimalizowania tego efektu można zastosować interpolację w węzłach Czebyszewa. W projekcie porównano jakość interpolacji z zastosowaniem węzłów Czebyszewa w stosunku do równoodległych węzłów. Węzły Czebyszewa w przedziale $[a, b]$ można obliczyć ze wzoru:
$ x_k = frac(a + b, 2) + frac(b - a, 2) cos(frac(2k - 1, 2k) pi), k = 1, ..., n $

== 2.2. Metoda interpolacji wykorzystująca funkcje sklejane trzeciego stopnia

Metoda funkcji sklejanych polega na wyznaczeniu wielomianu stałego stopnia na każdym z przedziałów, w których dzielony jest przedział interpolacji. Wielomiany są dobrane tak, aby były ciągłe oraz miały ciągłe pochodne pierwszego i drugiego rzędu. W projekcie zaimplementowano funkcje sklejane trzeciego stopnia, które są zdefiniowane na przedziale $[x_i, x_{i+1}]$ jako:
$ S_i (x) = a_i + b_i (x - x_i) + c_i (x - x_i)^2 + d_i (x - x_i)^3 $

#pagebreak()

Współczynniki $a_i, b_i, c_i, d_i$ można obliczyć na podstawie warunków brzegowych. W projekcie do implementacji metody funkcji sklejanych wykorzystano algorytm faktoryzacji LU, który pozwala na rozwiązanie układu równań liniowych.  Ponadto, metoda ta pozwala na interpolację funkcji o bardziej złożonym kształcie, co jest przydatne w przypadku analizy danych rzeczywistych.

$ cases(S_i (x_i) = y_i & "      i = 0, ..., n-1 wartość funkcji w węzłach",
S_i (x_(i+1)) = y_(i+1) & "      i = 0, ..., n-1 ciągłość funkcji w węzłach",
S''_0(x_0) = 0 & "      aerowanie drugiej pochodnej w punkcie początkowym",
S''_(n-1)(x_n) = 0 & "      zerowanie drugiej pochodnej w punkcie końcowym",
S'_i (x_i) = S'_(i-1)(x_i) & "      i = 1, ..., n-1 zerowanie pierwszej pochodnej w punktach",
S''_i (x_i) = S''_(i-1)(x_i) & "      i = 1, ..., n-1 zerowanie drugiej pochodnej w punktach") $

= 3. Dane wejściowe
W ramach projektu przeprowadzono testy na trzech różnych zestawach danych, które reprezentują profil wysokościowy terenu. Każdy z zestawów danych zawierał różną liczbę punktów oraz różne wartości wysokości. Dane  zostały pobrane z serwisu *`https://enauczanie.pg.edu.pl`*. Wybrane pliki użyte w projekcie z pliku *profile_wysokosciowe.zip* to:

- *`chelm.txt`* - plik reprezentujący profil wysokościowy terenu w okolicach i centrum Chełma w Gdańsku, punkty charakteryzują się małymi zmianami wysokości oraz małymi uskokami, widoczne sa płynne zmiany wysokości oraz niewielkie oscylacje miejscowe
- *`genoa_rapallo.txt`* - plik reprezentujący profil wysokościowy terenu miasta Genoa-Rapollo z dużym spadkiem pod koniec, punkty charakteryzują się dużymi zmianami wysokości, widoczne są duże oscylacje miejscowe
- *`tczew_starogard.txt`* - plik reprezentujący profil wysokościowy terenu w okolicach Tczewa i Starogardu, punkty charakteryzują stopniowym wzrostem wysokości, widoczne są duże oscylacje miejscowe

Wszystkie wykresy zawierają 512 punktów oraz zmienną liczbę punktów pomiarowych (6, 11, 15, 26, 52, 103). 

= 4. Analiza wyników metod interpolacji wielomianowej Lagrange’a i funkcji sklejanych trzeciego stopnia dla #underline[chelm.txt] - Trasa 1

#figure(
  image("plots/chelm/interpolation_6.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (6 punktów)]
)

#figure(
  image("plots/chelm/interpolation_11.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (11 punktów)]
)

#pagebreak()

Obserwując wykresy, szczególnie dla standardowej interpolacji wielomianowej Lagrange’a z równoodległymi węzłami, można zauważyć, że już dla 11 punktów pomiarowych interpolacja zaczyna stawać się mniej dokładna oraz zaczyna się pojawiać efekt Rungego na krańcach przedziału, który powoduje znaczne zaburzenia w interpolacji. W przypadku interpolacji wielomianowej Lagrange’a z węzłami Czebyszewa, efekt Rungego jest na tym samym poziomie o wiele mniej widoczny, co pozwala na uzyskanie lepszych wyników interpolacji. W przypadku funkcji sklejanych trzeciego stopnia, interpolacja jest bardzo dokładna nawet dla 11 punktów pomiarowych.

#figure(
  image("plots/chelm/interpolation_15.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (15 punktów)]
)

#figure(
  image("plots/chelm/interpolation_26.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (26 punktów)]
)


Przy 15 oraz 26 punktach pomiarowych, interpolacja jest bardzo nieprecyzyjna, a efekt Rungego uwydatnia się jeszcze bardziej i coraz bardziej się rozjeżdża. W przypadku funkcji z węzłami Czebyszewa, interpolacja jest mniej dokładna niż w przypadku funkcji sklejanych trzeciego stopnia, jednakże nadal pozwala na uzyskanie zadowalających wyników przy niewielkiej liczbie punktów pomiarowych. W przypadku interpolacji funkcjami sklejanymi trzeciego stopnia, interpolacja jest bardzo dokładna przy 26 punktach pomiarowych i praktycznie pokrywa się z oryginalnymi danymi z niewielkimi odchyleniami.

Zastosowanie 52 punktów pomiarowych pozwala na uzyskanie bardzo dokładnej interpolacji w przypadku wezłów Chebyshewa oraz funkcji sklejanych trzeciego stopnia. W przypadku interpolacji wielomianowej Lagrange’a z równoodległymi węzłami, efekt Rungego jest bardzo widoczny i każdy kolejny punkt pomiarowy powoduje coraz większe zaburzenia w interpolacji.

#figure(
  image("plots/chelm/interpolation_52.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (52 punkty)]
)

#figure(
  image("plots/chelm/interpolation_103.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (103 punkty)]
)

Obserwując kolejne wykresu zestawu danych *chelm.txt* można zauważyć, że liczba punktów pomiarowych ma znaczący wpływ na jakość interpolacji. W przypadku interpolacji wielomianowej Lagrange’a z równoodległymi węzłami, efekt Rungego jest bardzo widoczny już przy stosunkowo niewielkiej liczbie punktów pomiaru co ogranicza skuteczność tej metody oraz uniemożliwia poprawne odwzorowanie funkcji. Zastosowanie węzłów Czebyszewa pozwala na uniknięcie efektu Rungego, jednakże interpolacja jest mniej dokładna niż w przypadku funkcji sklejanych oraz w przypadku zbyt dużej ilości punktów pomiarowych nadal może istnieć prawdopodobnieństwo wystąpienia efektu Rungego. Najbardziej dokładna interpolacja została uzyskana przy użyciu funkcji sklejanych trzeciego stopnia.

Rozłożenie punktów pomiarowych ma znaczący wpływ na jakość interpolacji, jednakże w przypadku interpolacji funkcjami sklejanymi trzeciego stopnia, rozmieszczenie punktów pomiarowych nie ma większego znaczenia, a interpolacja zawsze jest dokładna. W przypadku zestawu danych Chełm, łagodność zmian wysokości pozwala na
uzyskanie dokładnej interpolacji nawet przy niewielkiej liczbie punktów pomiarowych kosztem pewnego uproszczenia profilu terenu.

= 5. Analiza wyników metod interpolacji wielomianowej Lagrange’a i funkcji sklejanych trzeciego stopnia dla #underline[genoa_rapallo.txt] - Trasa 2

#figure(
  image("plots/genoa_rapallo/interpolation_6.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (6 punktów)]
)

#figure(
  image("plots/genoa_rapallo/interpolation_11.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (11 punktów)]
)

#pagebreak()

Obserwując wykresy, szczególnie dla standardowej interpolacji wielomianowej Lagrange’a z równoodległymi węzłami, można zauważyć, że już dla pierwszego przypadku,czyli przy 6 punktach pomiarowych interpolacja zawodzi powodujac pojawienie się pierwszych oznak efektu Rungego. W przypadku 15 wezłów widoczne jest to jeszcze bardziej. Interpolacja z wezłami Czebyszewa pozwala na uzyskanie lepszych wyników, jednakże efekt Rungego jest nadal widoczny na krańcach przedziału przy zastosowaniu 15 wezłów. W przypadku funkcji sklejanych trzeciego stopnia, interpolacja na krańcu prawego przedziału przy 15 węzłach jest także zaburzona, jednakże jest to związane z gwałtownymi zmianami wysokości w tym obszarze.

#figure(
  image("plots/genoa_rapallo/interpolation_15.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (15 punktów)]
)

#figure(
  image("plots/genoa_rapallo/interpolation_26.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (26 punktów)]
)

Przy 26 punktach pomiarowych, interpolacja jest bardzo nieprecyzyjna, a efekt Rungego jest bardzo uwydatniony. Możemy zaobserwować duże oscylacje na krańcach przedziału, które powodują znaczne zaburzenia w interpolacji oraz występowanie harmoniczności. W przypadku funkcji z węzłami Czebyszewa, interpolacja o dziwo jest mniej dokładna niż przy zastosowaniu 15 punktów pomiarowych i końcowo wygląda gorzej niż interpolacja wielomianowa Lagrange’a z wezłami Czebyszewa dla 15 punktów pomiarowych. W przypadku interpolacji funkcjami sklejanymi trzeciego stopnia, interpolacja jest całkiem dokładna chociaż występują pewne ujemne zaburzenia na prawym krańcu przedziału.

#figure(
  image("plots/genoa_rapallo/interpolation_52.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (52 punkty)]
)

#figure(
  image("plots/genoa_rapallo/interpolation_103.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (103 punkty)]
)

Stosując 52 punkty pomiarowe, interpolacja jest bardzo dokładna w przypadku funkcji sklejanych trzeciego stopnia. W przypadku interpolacji wielomianowej Lagrange’a z równoodległymi węzłami, efekt Rungego stopniowo uwidacznia się jeszcze badziej. W przypadku interpolacji wielomianowej Lagrange’a z wezłami Czebyszewa, efekt Rungego zanika na lewym krańcu przedziału, jednakże na prawym krańcu przedziału jest nadal bardzo mocno widoczny. Przy zastosowaniu 103 punktów pomiarowych, interpolacja z węzłami Czebyszewa staje się nareszczie akceptowalna i zaczyna przypominać oryginalne dane. Jeśli chodzi o funkcje sklejane trzeciego stopnia, interpolacja jest bardzo dokładna i praktycznie pokrywa się z oryginalnymi danymi, chociaż występują pewne lekkie odchylenia w porównaniu do oryginalnych danych i poprzednich badanych danych wejściowych (*chelm.txt*).

Podsumowując analizę wyników interpolacji dla zestawu danych *genoa_rapallo.txt*, można zauważyć, że interpolacja wielomianowa Lagrange’a z równoodległymi węzłami jest bardzo podatna na efekt Rungego, który powoduje znaczne zaburzenia w interpolacji. Zastosowanie węzłów Czebyszewa pozwala na uniknięcie efektu Rungego, ale dopiero przy 103 punktach pomiarowych, jednakże interpolacja jest mniej dokładna niż w przypadku funkcji sklejanych trzeciego stopnia. Najbardziej dokładna interpolacja została uzyskana przy użyciu funkcji sklejanych trzeciego stopnia. Na podstawie tego zestawu danych można stwierdzić, że gwałtowne zmiany wysokości oraz duże oscylacje miejscowe mają znaczący wpływ na jakość interpolacji, a zastosowanie funkcji sklejanych trzeciego stopnia pozwala na uzyskanie najdokładniejszych wyników.
Widać także, że na przytoczonym wykresie interpolacja wielomianowa Lagrange’a z równoodległymi węzłami nie radzi sobie z gwałtownymi zmianami wysokości, dając niezadowalające wyniki interpolacji. Czebyszewa również daje niezadowalające wyniki. Widać oscylację na końcu profilu terenu, z uwagi na duża liczbę gwałtownych zmian wysokości. 

= 6. Analiza wyników metod interpolacji wielomianowej Lagrange’a i funkcji sklejanych trzeciego stopnia dla #underline[tczew_starogard.txt] - Trasa 3

#figure(
  image("plots/tczew_starogard/interpolation_6.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (6 punktów)]
)

#pagebreak()

#figure(
  image("plots/tczew_starogard/interpolation_11.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (11 punktów)]
)

Obserwując wykresy, szczególnie dla standardowej interpolacji wielomianowej Lagrange’a z równoodległymi węzłami, można zauważyć, że już dla 6 punktów pomiarowych interpolacja staje się niedokładna na lewym krańcu przedziału. W przypadku 11 punktów pomiarowych, interpolacja jest bardziej dokładna, jednakże efekt Rungego zaczyna się uwydatniać na obu krańcach. W przypadku pozostałych dwóch metod interpolacji na 11 wezłach, interpolacja jest zadowalająca i pozwala na uzyskanie całkiem dokładnych wyników. 

#figure(
  image("plots/tczew_starogard/interpolation_15.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (15 punktów)]
)

#figure(
  image("plots/tczew_starogard/interpolation_26.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (26 punktów)]
)

Przy 26 punktach pomiarowych, interpolacja jest bardzo nieprecyzyjna, a efekt Rungego jest bardzo uwydatniony na prawym i lewym krańcu. Możemy zaobserwować duże oscylacje na krańcach przedziału, które powodują znaczne zaburzenia w interpolacji oraz występowanie harmoniczności. W przypadku funkcji z węzłami Czebyszewa, interpolacja zaczyna prawidłowo pokrywać się z oryginalnymi danymi, jednakże nie jest to jeszcze idealna interpolacja. W przypadku interpolacji funkcjami sklejanymi trzeciego stopnia, interpolacja jest bardzo dokładna i praktycznie pokrywa się z oryginalnymi danymi z niewielkimi odchyleniami. Najgorsze wyniki uzyskano dla interpolacji wielomianowej Lagrange’a z równoodległymi węzłami, gdzie efekt Rungego jest bardzo widoczny i każdy kolejny punkt pomiarowy powoduje coraz większe zaburzenia w interpolacji. Najbardziej jest to widoczne na najniższym etapie wzniesienia dla banego profilu terenu.

#figure(
  image("plots/tczew_starogard/interpolation_52.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (52 punkty)]
)

#figure(
  image("plots/tczew_starogard/interpolation_103.png"),
  caption: [Interpolacja wielomianowa Lagrange’a i funkcji sklejanych trzeciego stopnia (103 punkty)]
)

Stosując 52 punkty pomiarowe, interpolacja jest bardzo dokładna w przypadku funkcji sklejanych trzeciego stopnia. W przypadku interpolacji wielomianowej Lagrange’a z równoodległymi węzłami, efekt Rungego stopniowo uwidacznia się jeszcze badziej. W przypadku interpolacji wielomianowej Lagrange’a z wezłami Czebyszewa zaczynamy uzyskiwać jeszcze bardziej dokładne wyniki. Przy zastosowaniu 103 punktów pomiarowych, interpolacja z węzłami Czebyszewa staje się bardzo dobra, a funkcja sklejana trzeciego stopnia jest bardzo dokładna i praktycznie pokrywa się z oryginalnymi danymi. Interpolacja wielomianowa Lagrange’a z równoodległymi węzłami jest ponownie bardzo podatna na efekt Rungego, który powoduje znaczne zaburzenia w interpolacji. 

Podsumowując analizę wyników interpolacji dla zestawu danych *tczew_starogard.txt*, można zauważyć, że interpolacja wielomianowa Lagrange’a z równoodległymi węzłami jest bardzo podatna na efekt Rungego, który powoduje znaczne zaburzenia w interpolacji. Zastosowanie węzłów Czebyszewa pozwala na uniknięcie efektu Rungego, jednakże interpolacja jest mniej dokładna niż w przypadku funkcji sklejanych trzeciego stopnia. Najbardziej dokładna interpolacja została uzyskana przy użyciu funkcji sklejanych trzeciego stopnia. Na podstawie tego zestawu danych można stwierdzić, że metody interpolacji są wrażliwe na gwałtowne zmiany wysokości oraz duże oscylacje miejscowe, które mają znaczący wpływ na jakość interpolacji. Najgorzej sobie poradziły z początkiem wzniesienia profilu terenu.

= 7. Porównanie metod interpolacji i podsumowanie

Podsumowując wyniki przeprowadzonych testów, można stwierdzić, że metoda interpolacji wielomianowej Lagrange’a z równoodległymi węzłami jest bardzo podatna na efekt Rungego, który powoduje znaczne zaburzenia w interpolacji. Zastosowanie węzłów Czebyszewa pozwala na uniknięcie efektu Rungego, jednakże interpolacja jest mniej dokładna niż w przypadku funkcji sklejanych trzeciego stopnia. Najbardziej dokładna interpolacja została uzyskana przy użyciu funkcji sklejanych trzeciego stopnia. 

W przypadku wszystkich zestawów danych można zauważyć, że interpolacja wielomianowa Lagrange’a z równoodległymi węzłami jest najgorszym wyborem jeśli chcemy uzyskać dokładne wyniki interpolacji przy aproksymacji profilu wysokościowego terenu. Zastosowanie węzłów Czebyszewa pozwala na uniknięcie efektu Rungego w niektórych przypadkach i uzyskanie lepszych wyników, jednakże nie jest ona bezbłędna. Najbardziej dokładne wyniki uzyskano przy użyciu funkcji sklejanych trzeciego stopnia, które pozwoliły na uzyskanie dokładnej interpolacji nawet przy niewielkiej liczbie punktów pomiarowych.

Interpolacja z wezłami Czebyszewa wymaga dodatkowo wyboru nierównomiernie rozmieszczonych węzłów, co jest trudne do osiągnięcia w praktyce. W plikach z danymi wejściowymi znajdowało się aż 512 punktów, co pozwalało na odpowiednie dobranie węzłów Czebyszewa. W przypadku mniejszej ilości danych, nie byłoby to możliwe. Z tego
powodu metoda Lagrange’a nie jest odpowiednia do interpolacji profilu wysokościowego.

Metoda funkcji sklejanych pozwala na uzyskanie bardzo dokładnych wyników interpolacji. Jest ona odporna na efekt Rungego, co pozwala na zwiększanie ilości węzłów bez obaw o pogorszenie wyników. Metoda funkcji sklejanych spełnia swoje zadanie dla danych równomiernie rozmieszczonych, co jest idealne w przypadku profilu wysokościowego.

Podsumowując, pomimo trudniejszej implementacji, metoda funkcji sklejanych trzeciego stopnia jest najlepszym wyborem do interpolacji profilu wysokościowego terenu. Pozwala na uzyskanie dokładnych wyników interpolacji nawet przy niewielkiej liczbie punktów pomiarowych oraz jest odporna na efekt Rungego. Metoda ta pozwala na interpolację funkcji o bardziej złożonym kształcie, co jest przydatne w przypadku analizy danych rzeczywistych. Dodatkowo mozliwe jest zastosowanie tej metody do interpolacji danych o różnym rozkładzie punktów pomiarowych bez utraty dokładności interpolacji przy zwiększaniu ilości punktów pomiarowych.