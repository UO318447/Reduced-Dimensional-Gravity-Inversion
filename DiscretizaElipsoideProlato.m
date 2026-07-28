%*******************************************************************************
% Function:  [puntos,dvr,t] = DiscretizaElipsoideProlato(a,b,pe,tmax)
%
% Propósito: Discretiza un elipsoide prolato en un conjunto de cubos
%
% Entradas:  - a: Semieje mayor del elipsoide, coincidente con el eje Y
%            - b: Semieje menor del elipsoide, coincidente con los ejes X y Z
%                 (la elipse gira en torno al eje Y para generar el elipsoide)
%            - pe: Porcentaje máximo de error en el volumen de la aproximación
%                  del elipsoide discretizado, en tanto por CIENTO
%            - tmax: Tiempo máximo, en segundos, de ejecución de la función.
%                    Superado este tiempo la función termina aunque no se haya
%                    alcanzado el nivel de error en volumen solicitado. Si se
%                    pasa el valor 0 la función continúa con el criterio del
%                    volumen hasta que se alcance la tolerancia establecida
%
% Salidas:   - puntos: Matriz de tres columnas con las coordenadas de los
%                      centros de los cubos en los que se discretiza el
%                      elipsoide. Cada fila es un punto
%                      - Col. 1: Coordenada X del centro del cubo
%                      - Col. 2: Coordenada Y del centro del cubo
%                      - Col. 3: Coordenada Z del centro del cubo
%            - dvr: Diferencia relativa, en tanto por CIENTO) entre el volumen
%                   del elipsoide discretizado y el real
%            - t: Tiempo transcurrido de ejecución, en segundos
%
% Nota 1: No se realiza ninguna comprobación sobre los argumentos de entrada
% Nota 2: Para elipsoides con excentricidades muy grandes, por ejemplo, el
%         semieje mayor 100 veces más grande que el menor, la función puede
%         ralentizarse mucho y generar un número de cubos excesivo
%
% Historia:  05-07-2019: Creación de la función
%                        José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************

function [puntos,dvr,t] = DiscretizaElipsoideProlato(a,b,pe,tmax)

%Marca inicial de tiempo
tic;
%Comprobamos degeneraciones
if (a<=0.0)||(b<=0.0)||(a<b)
    error('Uno o ambos semiejes del elipsoide es o son <= 0.0, o a<b');
end
%Volumen del elipsoide
ve = 4.0/3.0*pi*a*b^2;
%Cubo concéntrico con el elipsoide y de arista igual al semieje menor
lado = b;
%Entramos en un bucle infinito
while 1
    %Coordenadas de los centros de los cubos que rellenan el paralelepípedo en
    %el que se inscribe el elipsoide
    cx = [fliplr(-lado:-lado:-b) 0.0:lado:b];
    cy = [fliplr(-lado:-lado:-a) 0.0:lado:a];
    cz = cx;
    %Número de cubos en cada eje
    ncx = length(cx);
    ncy = length(cy);
    ncz = length(cz);
    %Variable de salida
    puntos = zeros(ncx*ncy*ncz,3);
    %Genero las coordenadas de los centros de los cubos
    [x,y,z] = meshgrid(cx,cy,cz);
    puntos(:,1) = x(:);
    puntos(:,2) = y(:);
    puntos(:,3) = z(:);
    %Calculamos la ecuación del elipsoide para todos los puntos
    ecelip = (puntos(:,2).^2)./(a^2)+(puntos(:,1).^2+puntos(:,3).^2)./(b^2);
    %Buscamos los puntos dentro o en la superficie
    pds = ecelip<=1.0;
    %Calculamos el volumen compuesto por los cubos seleccionados
    vc = sum(pds)*lado^3;
    %Ciferencia relativa entre el volumen real y el calculado
    dvr = abs(ve-vc)/ve*100.0;
    %Comprobamos si la diferencia es menor que la tolerancia establecida
    t = toc;
    if (dvr<=pe)||((tmax~=0.0)&&(t>tmax))
        %Nos quedamos con los puntos del elipsoide
        puntos = puntos(pds,:);
        %Salimos de la función
        break;
    else
        %Hacemos que en el eje más corto quepa un cubo más
        lado = (cx(end)-cx(1))/ncx;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright (c) 2019-2024, J.L.G. Pallero, jgpallero@gmail.com
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
