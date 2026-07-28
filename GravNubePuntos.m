%*******************************************************************************
% Function:  [g] = GravNubePuntos(atractores,atraidos)
%
% Propósito: Calcula la atracción gravitacional de una nube de puntos sobre un
%            conjunto de puntos
%
% Entradas:  - atractores: Matriz de tres o cuatro columnas con la nube de
%                          puntos atractores. Cada fila es un punto:
%                          - Col. 1: Coordenada X del punto, en metros
%                          - Col. 2: Coordenada Y del punto, en metros
%                          - Col. 3: Coordenada Z del punto, en metros
%                          - Col. 4: Masa del punto, en kg
%                          Si la matriz tiene tres columnas se asume una masa
%                          igual a la unidad para todos los puntos
%            - atraidos: Matriz de tres columnas que contiene el conjunto de
%                        puntos atraídos. Cada fila es un punto:
%                        - Col. 1: Coordenada X del punto, en metros
%                        - Col. 2: Coordenada Y del punto, en metros
%                        - Col. 3: Coordenada Z del punto, en metros
%
% Salidas:   - g: Matriz de tres columnas con un número de elementos igual al
%                 número de filas de la matriz 'atraidos'. Cada fila contiene
%                 las componentes del vector atracción gravitacional del
%                 conjunto de puntos atractores sobre el punto atraído
%                 correspondiente. Las columnas son:
%                 - Col. 1: Componente X del vector atracción, en m/s^2
%                 - Col. 2: Componente Y del vector atracción, en m/s^2
%                 - Col. 3: Componente Z del vector atracción, en m/s^2
%
% Nota: No se realiza ninguna comprobación sobre los argumentos de entrada
%
% Historia:  28-09-2024: Creación de la función
%                        José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************

function [g] = GravNubePuntos(atractores,atraidos)

%Constante de gravitación universal
G = 6.67430e-11;
%Dimensiones de las numes de puntos
[npa,col] = size(atractores);
npA = size(atraidos,1);
%Si no hay masas, las añadimos
if col==3
    atractores = [atractores ones(npa,1)];
end
%Variable de salida
g = zeros(npA,3);
%Variable auxiliar
v = zeros(npa,3);
%Recorremos los puntos atraídos
for i=1:npA
    %Vectores de la masa atraída a los puntos atractores
    v(:,1) = atraidos(i,1)-atractores(:,1);
    v(:,2) = atraidos(i,2)-atractores(:,2);
    v(:,3) = atraidos(i,3)-atractores(:,3);
    %Módulo al cuadrado
    d2 = v.^2;
    d2 = sum(d2')';
    %Módulo
    d = sqrt(d2);
    %Parte común de la fórmula
    d = -G*atractores(:,4)./d;
    %Componentes de la atracción
    g(i,1) = g(i,1)+sum(d.*v(:,1)./d2);
    g(i,2) = g(i,2)+sum(d.*v(:,2)./d2);
    g(i,3) = g(i,3)+sum(d.*v(:,3)./d2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 2024-2025, J.L.G. Pallero, jgpallero@gmail.com
%
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%- Redistributions of source code must retain the above copyright notice, this
%  list of conditions and the following disclaimer.
%- Redistributions in binary form must reproduce the above copyright notice,
%  this list of conditions and the following disclaimer in the documentation
%  and/or other materials provided with the distribution.
%- Neither the name of the copyright holders nor the names of its contributors
%  may be used to endorse or promote products derived from this software without
%  specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
%INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
%OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
%ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
