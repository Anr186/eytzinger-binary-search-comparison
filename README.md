# Сравнительный анализ эффективности алгоритмов поиска в непрерывно размещенных в памяти массивах  Терентьев А.В., Коротков К.К. 

*Исследование выполнено в Сибирском государственном университете телекоммуникаций и информатики*

## 📌 Авторы
- **Терентьев Андрей Викторович** (andrey.terent1ev@mail.ru)
- **Коротков Кирилл Константинович** (karton1228@gmail.com)  
*Научный руководитель: д.т.н. профессор Курносов М.Г.*

---

## 🔍 Цель исследования
Сравнение производительности алгоритмов поиска:
- **Классический двоичный поиск** (BS/BL/Prefetch)
- **Поиск Eytzinger** (Eyt/EytBL)  
на архитектурах **x86-64 (Xeon Gold)** и **RISC-V (SpacemiT K1)**.

---

## 📊 Ключевые результаты
### 1. Эффективность алгоритмов
| Алгоритм       | Ускорение (L1/L2) | Преимущества                          |
|----------------|-------------------|---------------------------------------|
| Branchless (BL)| 1.5–1.75×         | Устранение условных переходов         |
| Eytzinger      | 1.75–2×           | Лучшая локальность данных в кэше L3   |

### Графики времени выполнения
![] (img/1.png)
![] (img/2.png)
![] (img/3.png)
![] (img/4.png)

*Зависимость времени поиска от размера массива*
### 2. Влияние кэш-памяти
- **L1/L2**: Branchless-оптимизации дают максимальный эффект.
- **L3+**: Eytzinger-размещение сокращает промахи кэша.  
- **RISC-V**: Отсутствие L3 кэша снижает преимущество Eytzinger.

---

## 🛠 Технические детали
### Реализованные алгоритмы
```cpp
// Пример branchless-поиска
IdxType branchfree_search(ElemType x) const {
    const ElemType *base = a; 
    IdxType n = this->n;
    while (n > 1) {
        const IdxType half = n / 2;
        base = (base[half] < x) ? &base[half] : base;
        n -= half;
    }
    return (*base < x) + base - a;
}
// Пример branchless-поиска с явной предвыборкой
IdxType branchfree_pref_search(ElemType x) const {
    const ElemType *base = a; IdxType n = this->n;
    while (n > 1) {
      const IdxType half = n / 2;
      builtin_prefetch(&base[half/2]);
      builtin_prefetch(&base[half + half/2]);
      base = (base[half] < x)? &base[half]: base;
      n -= half;
    }
    return (*base < x) + base - a;
}

// Пример branchless Eytzinger-поиска
IdxType eytz_branchless_search(ElemType x) const {
    IdxType i = 0;
    while (i < n) {
      i = (x <= a[i]) ? (2*i + 1) : (2*i + 2);
    }    
    IdxType j = (i+1) >> __builtin_ffs(~(i+1));
    return (j == 0) ? n : j-1;
}
```
Доклад на всероссийской научно-технической конференции "Обработка информации и математическое моделирование (ОИиММ-2025)", Новосибирск, 2025.
