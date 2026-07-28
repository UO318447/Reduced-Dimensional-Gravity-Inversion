%*******************************************************************************
% Function:  [g] = GravElipsoideProlatoY(elipsoide,puntos)
%
% Propósito: Calcula la atracción gravitacional de un elipsoide prolato
%            orientado según el eje Y sobre un conjunto de puntos exteriores
%
% Entradas:  - elipsoide: Vector de dos o tres elementos con la definición del
%                         elipsoide de trabajo, que estará centrado en el origen
%                         de coordenadas
%                         - Pos. 1: Semieje mayor del elipsoide, en metros, que
%                                   estará orientado en la dirección del eje Y
%                         - Pos. 2: Semieje menor del elipsoide, en metros, que
%                                   corresponderá a las direcciones X y Z
%                         - Pos. 3: Densidad del elipsoide, en kg/m^3. Si se
%                                   omite este elemento se asumirá que la
%                                   densidad es igual a la unidad
%            - puntos: Matriz de tres columnas que contiene el conjunto de
%                      puntos atraídos, referidos al origen de coordenadas
%                      centrado en el elipsoide. Cada fila es un punto:
%                      - Col. 1: Coordenada X del punto, en metros
%                      - Col. 2: Coordenada Y del punto, en metros
%                      - Col. 3: Coordenada Z del punto, en metros
%
% Salidas:   - g: Matriz de tres columnas con un número de elementos igual al
%                 número de filas de la matriz 'puntos'. Cada fila contiene las
%                 componentes del vector atracción gravitacional del elipsoide
%                 sobre el punto correspondiente. Las columnas son:
%                 - Col. 1: Componente X del vector atracción, en m/s^2
%                 - Col. 2: Componente Y del vector atracción, en m/s^2
%                 - Col. 3: Componente Z del vector atracción, en m/s^2
%
% Nota 1: No se realiza ninguna comprobación sobre los argumentos de entrada
% Nota 2: Formulación tomada de W.D. MacMillan (1958), The Theory of the
%         Potential. Dover Publications, 1a edición. ISBN: 978-0486604862
%
% Historia:  28-09-2024: Creación de la función
%                        José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************

function [g] = GravElipsoideProlatoY(elipsoide,puntos)

%Constante de gravitación universal
G = 6.67430e-11;
%Semiejes del elipsoide y densidad
a = elipsoide(1);
c = elipsoide(2);
ro = 1.0;
if length(elipsoide)>2
    ro = elipsoide(3);
end
%Comprobamos degeneraciones
if (a<=0.0)||(c<=0.0)||(a<c)
    error('Uno o ambos semiejes del elipsoide es o son <= 0.0, o a<c');
end
%Potencias de las coordenadas de los puntos de trabajo
x2 = puntos(:,1).^2;
y2 = puntos(:,2).^2;
z2 = puntos(:,3).^2;
%Potencias de los semiejes del elipsoide
a2 = a^2;
c2 = c^2;
%Parámetros kappa
B = a2+c2-x2-y2-z2;
C = a2*(c2-x2-z2)-c2*y2;
k = (-B+sqrt(B.^2-4.0*C))./2.0;
%Parte común a todas las componentes de la atracción
A = 2.0*pi*G*ro*c2*a/(a2-c2);
B = asinh(sqrt((a2-c2)./(c2+k)))./sqrt(a2-c2);
C = sqrt(a2+k)./(c2+k);
%Variable de salida
g = zeros(size(puntos,1),3);
%Componente X
g(:,1) = A*puntos(:,1).*(B-C);
%Componente Y
g(:,2) = A*puntos(:,2).*(-2.0*B+2.0./sqrt(a2+k));
%Componente Z
g(:,3) = A*puntos(:,3).*(B-C);
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
