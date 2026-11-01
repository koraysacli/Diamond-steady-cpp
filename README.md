# Diamond: Mach 3, geniş alan, steady OpenFOAM

OpenFOAM **Foundation 14**, `shockFluid`, 2D inviscid (Euler) akış. Hücum açısı **0°**, irtifa **80.000 ft**. Profil kord uzunluğu 0,2 m, tam kalınlığı 0,02 m; tek hücrelik düzlem dışı kalınlık 0,001 m. Profil, diamond.stl dosyasının dört köşeli kesitiyle aynıdır; ağ doğrudan bu kesite oturtulur.

Akış alanı **x = −1…4 m, y = −2,5…2,5 m**. İlk 24.000 hücreli vakanın 0,6 × 0,4 m alanının yaklaşık 104 katı. Profil etrafında sıklaşan **1.200.000 altı yüzlü hücre**, **20 MPI süreci**. Akış alanı büyüdüğü için bütün alanı gösterirken profil ekranda küçük görünür.

## Kod nerede?

```text
src/CaseConfig.hpp      ← Mach, irtifa, açı, boyutlar; atmosfer if/else formülleri
src/generateCase.cpp    ← for/if/else ile ağ, sınır koşulu ve sözlük üretimi
src/inletVelocity.C     ← OpenFOAM içinde derlenen C++ hız koşulu
src/diagnostics.C       ← C++ basınç kuvveti, Cl/Cd ve kütle dengesi
src/shockReference.cpp  ← döngüyle çözülen eğik şok ve genleşme referansı
0/, system/, constant/  ← C++ programının ürettiği OpenFOAM giriş dosyaları
monitor.py              ← kayıtları okuyan Python izleyici (çözücüye dokunmaz)
Allrun                  ← derleme, ağ, kontrol, paralel çözüm, reconstruction
```

Uzun hız tabloları yoktur. OpenFOAM kendi sözlük dosyalarını gerektirir; bunları C++ üretir. Hız koşulu `codedFixedValue` ile çalışma sırasında ayrıca derlenir. C++ kaynaklarını değiştirmek çalışan analize kendiliğinden uygulanmaz. Değişiklikten sonra **yeni bir vaka dizini** üretin; mevcut mesh/sonuçların üzerine yazma koruması vardır.

```bash
mkdir -p build
g++ -std=c++17 -O2 -Wall -Wextra "$PWD/src/generateCase.cpp" -o build/generateCase
build/generateCase /yeni/vaka/klasoru
```

## Çalıştırma ve görüntüleme

Ubuntu'da Foundation 14 ve `g++` kurulu olmalı. Depo ana dizininde:

```bash
bash Allrun > log.Allrun 2>&1
```

Ayrı bir terminalde, aynı dizinde:

```bash
watch -n 10 python3 monitor.py
```

Bu ekranda Ctrl+C sadece izleyiciyi kapatır. **Allrun çalışan terminalde Ctrl+C çözümü keser.** Kaydederek durdurmak için OpenFOAM ortamı yüklü ayrı terminalde:

```bash
foamDictionary system/controlDict -entry stopAt -set writeNow
```

Devam etmeden önce `stopAt` değerini `endTime` yapın; `startFrom latestTime` son kayıttan başlar. Analizi iki kez aynı anda başlatmayı `flock` engeller.

```bash
paraview --no-mpi diamond.foam
```

Çözüm sürerken ParaView'da **Case Type → Decomposed Case**, `internalMesh`, Apply ve son kayıtlı zamanı seçin. Tamamlanınca `reconstructPar -latestTime` çalışır; **Reconstructed Case** ile açılabilir. Bu vaka başka klasördeki eski `diamond.foam` dosyasından ayrıdır. `p`, `U` ve `T` alanları mevcuttur. `Time` değeri aşağıdaki LTS açıklamasına göre okunmalıdır.

## Steady ne anlama geliyor?

`localEuler` yerel zaman adımıyla kararlı çözüme ilerler. Ekrandaki `Time = 1000s` yazısı burada **1000 sözde-zaman iterasyonu** demektir; 1000 saniyelik fiziksel transient değildir. Sabit sınır koşulları vardır. Maksimum 20.000 iterasyon otomatik bir yakınsama garantisi değildir.

Kurganov akı şeması, `Minmod/MinmodV` sınırlayıcı ve yerel `maxCo = 0.2` kullanılır. İlk `vanLeer`, Co=0.5 denemesi büyüyen salınım ve negatif sıcaklık kaynaklı sayısal hata verdi; o deneme sonuç kabul edilmedi. Son tercih daha sönümlüdür; şok kalınlığı ve sürükleme için ağ/şema duyarlılığı yine incelenmelidir.

Doğrusal çözücüde `diagonal` için görülen sıfır residual **yakınsama kanıtı değildir**. `monitor.py`, en az 5000 iterasyon ve son 2000 iterasyonluk pencereyi denetler: Cd değişimi <1e−5, |Cl|<1e−5, normalize net kütle akısı <1e−4, 20 adım arası hacim ağırlıklı ΔCp RMS <1e−5, ön/arka yüzey ortalama basınç oranlarının değişimi <0,001. Bunlar mühendislik eşikleri; alan çok büyük olduğu için hacim ağırlıklı değişim yerel hatayı gizleyebilir. Şok konumu, yüzey Cp ve ağ bağımsızlığı ayrıca incelenmelidir.

## Kuvvetler ve şok

`F = Σ (p − p∞) Sf`, burada `Sf` akışkandan katıya yönelir. Referans alanı **kord × düzlem dışı kalınlık = 0,0002 m²**, dinamik basınç `q∞ = ρ∞ U∞²/2`. Akış yönündeki kuvvetten Cd, dik yönden Cl hesaplanır. Her MPI sürecinin katkısı toplanır.

Modelde **viskozite sıfır**: Cd yalnızca basınç/dalga sürüklemesidir; gerçek toplam sürükleme veya sınır tabaka ayrılması olarak yorumlanmaz. 0° ve simetrik profil için Cl sıfıra yaklaşmalıdır. Gerçek viskoz çözüm için duvar ağı, Reynolds sayısı ve uygun türbülans/laminer model ayrı bir çalışma gerektirir.

`validation/theory.txt` analitik referanstır, CFD sonucu değildir. Zayıf eğik şok + Prandtl–Meyer genleşmesi: β≈23,70098°, Cd≈0,01429478, Cl=0, ön yüzey p/p∞≈1,529866, arka yüzey≈0,629295. Şok bağıntıları için [NASA eğik şok açıklaması](https://www.grc.nasa.gov/WWW/BGH/oblique.html). Çözümden ayrıca şok açısı, basınç sıçraması, Mach alanı ve yüzey Cp kontrol edilmelidir.

## Doğrulama ve paylaşım

`checkMesh`: 1.200.000 hücre, pozitif hacimler, maksimum non-orthogonality 5,71°, skewness 0,244, **Mesh OK**. 12.000 hücrelik ayrı testte C++ sınır koşulu/izleyici derlendi; 40 adımda seri ve 2 süreçli Cd/Cl farkı 1e−6 altında, temiz kapanış doğrulandı. Bu kısa test yakınsama veya ağ bağımsızlığı testi değildir.

GitHub'a kaynaklar, küçük giriş dosyaları ve doğrulama özeti gönderilir. `.gitignore`: mesh, processor dizinleri, zaman sonuçları, loglar, dynamicCode ve derlenmiş ürünleri dışarıda bırakır. Üretilmiş sonuçlar yerelde korunur.
