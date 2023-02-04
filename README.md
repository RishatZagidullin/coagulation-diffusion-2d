## Code description

Код решает уравнение 2Д диффузии и коагуяции. Для компилляции достаточно иметь gcc. Есть мейкфайл, в конце должен выдать файл исполнения под названием **tet.exe**. Картинки (в формате .ppm) после исполнения кода будут храниться в **run_folder/imgs**.

Аргументы командной строки: ``./tet.exe 80 50 0`` - ``80`` - x-координата источника мономеров, ``50`` - y-координата источника мономеров, ``0`` - включить ли расчет коагуляции (0 - нет, 1 - да).

На двумерной плоскости задается в любой точке источник мономеров, и дальше делается расчет.

Что можно менять:
* граничные условия;
* добавить перенос;
* расчетная область;
* размеры картинки;
* ядро коагуляции;
* коэффициенты в уравнении;
* координаты источника;
* шаги по времени/пространству/размеру частиц.
