%*******************************************************************************
% Function:  [g,esup,dvr,edisc] = GravElipProlato(elip,puntos,algoritmo,pe,tmax)
%
% Propósito: Calcula la atracción gravitacional de un conjunto de elipsoides
%            prolatos sobre un conjunto de puntos exteriores a todos ellos
%
% Entradas:  - elip: Matriz de ocho columnas donde cada fila es la definición de
%                    un elipsoide
%                    - Col. 1: Coordenada X del centro del elipsoide, en metros
%                    - Col. 2: Coordenada Y del centro del elipsoide, en metros
%                    - Col. 3: Coordenada Z del centro del elipsoide, en metros
%                    - Col. 4: Semieje mayor del elipsoide, en metros
%                    - Col. 5: Semieje menor del elipsoide, en metros
%                    - Col. 6: Acimut del semieje mayor del elipsoide, medido
%                              como en topografía, es decir, positivo en sentido
%                              horario y con origen en el eje Y (asimilado al
%                              norte). En radianes
%                    - Col. 7: Inclinación del semieje mayor del elipsoide,
%                              considerado positivo hacia arriba y negativo
%                              hacia abajo. En radianes
%                    - Col. 8: Densidad del elipsoide, en kg/m^3
%            - puntos: Matriz de tres columnas que contiene el conjunto de
%                      puntos atraídos. Cada fila es un punto:
%                      - Col. 1: Coordenada X del punto, en metros
%                      - Col. 2: Coordenada Y del punto, en metros
%                      - Col. 3: Coordenada Z del punto, en metros
%            - algoritmo: Identificador del algoritmo a emplear para los
%                         cálculos. Dos posibles valores:
%                         - 0: Se utiliza la formulación exacta de MacMillan
%                         - Distinto de 0: Se utiliza la discretización de los
%                              elipsoides en cubos, que se tratarán como puntos
%                              en el cálculo de la atracción
%                         Por omisión, este argumento vale 0
%            - pe: Porcentaje máximo de error en el volumen de la aproximación
%                  del elipsoide discretizado, en tanto por CIENTO. Este
%                  argumento sólo se tiene en cuenta si el valor de 'algoritmo'
%                  es distinto de 0
%            - tmax: Tiempo máximo, en segundos, de ejecución de la función.
%                    Superado este tiempo la función termina aunque no se haya
%                    alcanzado el nivel de error en volumen solicitado. Si se
%                    pasa el valor 0 la función continúa con el criterio del
%                    volumen hasta que se alcance la tolerancia establecida
%
% Salidas:   - g: Matriz de tres columnas con un número de elementos igual al
%                 número de filas de la matriz 'puntos'. Cada fila contiene las
%                 componentes del vector atracción gravitacional del conjunto de
%                 elipsoides sobre el punto correspondiente. Las columnas son:
%                 - Col. 1: Componente X del vector atracción, en m/s^2
%                 - Col. 2: Componente Y del vector atracción, en m/s^2
%                 - Col. 3: Componente Z del vector atracción, en m/s^2
%            - esup: Cell array con las coordenadas de una serie de puntos en la
%                    superficie de los elipsoides. Cada elemento refiere al
%                    elipsoide de la fila correspondiente del argumento 'elip' y
%                    es una matriz de tres columnas con las coordenadas de los
%                    centros de los cubos en los que se ha discretizado. Cada
%                    fila es un punto y las columnas son:
%                     - Col. 1: Coordenada X del punto, en metros
%                     - Col. 2: Coordenada Y del punto, en metros
%                     - Col. 3: Coordenada Z del punto, en metros
%                     Las coordenadas anteriores se refieren al sistema GLOBAL,
%                     es decir, al que vienen referidos los puntos atraídos
%            - dvr: Vector con la diferencia relativa, en tanto por CIENTO,
%                   entre el volumen del elipsoide correspondiente a su posición
%                   discretizado y el real. Este argumento de salida no tiene
%                   sentido si 'algoritmo' es igual a 0, por lo que en tal caso
%                   se devuelve un vector vacío
%            - edisc: Cell array con la discretización de los elipsoides. Cada
%                     elemento refiere al elipsoide de la fila correspondiente
%                     del argumento 'elip' y es una matriz de tres columnas con
%                     las coordenadas de los centros de los cubos en los que se
%                     ha discretizado. Cada fila es un punto y las columnas son:
%                     - Col. 1: Coordenada X del centro del cubo, en metros
%                     - Col. 2: Coordenada Y del centro del cubo, en metros
%                     - Col. 3: Coordenada Z del centro del cubo, en metros
%                     Las coordenadas anteriores se refieren al sistema GLOBAL,
%                     es decir, al que vienen referidos los puntos atraídos
%                     Este argumento de salida no tiene sentido si 'algoritmo'
%                     esigual a 0, por lo que en tal caso se devuelve un cell
%                     array vacío
%
% Nota 1: La formulación para el cálculo exacto ha sido tomada de
%         W.D. MacMillan (1958), The Theory of the Potential
%         Dover Publications, 1a edición. ISBN: 978-0486604862
% Nota 2: El convenio de signos es el que se desprende de considerar las
%         componentes del vector intensidad como Ix=dV/dx, Iy=dV/dy e Iz=dV/dz,
%         donde V=Gm/r, donde V es el potencial gravitacional, m la masa
%         atractora y r la distancia que la separa del punto de trabajo
%
% Historia:  28-09-2024: Creación de la función
%                        José Luis García Pallero, jgpallero@gmail.com
%*******************************************************************************

function [g,esup,dvr,edisc] = GravElipProlato(elip,puntos,algoritmo,pe,tmax)

%Número de argumentos de entrada
if (nargin<2)||((nargin>=3)&&(nargin<5))
    error('Número incorrecto de argumentos de entrada');
elseif nargin==2
    algoritmo = 0;
    pe = 1.0;
    tmax = 0.0;
end
%Dimensiones de las matrices de elipsoides y puntos
[nElip,colElip] = size(elip);
[nP,colP] = size(puntos);
if (colElip<8)||(colP<3)
    error(['Las dimensiones de las matrices ''elip'' y ''puntos'' son ',...
           'incorrectas']);
end
%Valores de los dos últimos argumentos
if pe<=0.0
    error('El argumento ''pe'' tiene un valor incorrecto');
end
if tmax<0.0
    error('El argumento ''tmax'' tiene un valor incorrecto');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variables de salida
g = zeros(nP,3);
esup = {};
dvr = [];
edisc = {};
if nargout>1
    esup = cell(nElip,1);
    if (nargout>2)&&(algoritmo~=0)
        dvr = zeros(nElip,1);
        if (nargout>3)&&(algoritmo~=0)
            edisc = cell(nElip,1);
        end
    end
end
%Recorro el número de elipsoides
for i=1:nElip
    %Coordenadas del centro del elipsoide
    xe = elip(i,1);
    ye = elip(i,2);
    ze = elip(i,3);
    %Ángulos de rotación
    acim = elip(i,6);
    inc = elip(i,7);
    %Matriz de rotación para llevar los puntos al sistema del elipsoide
    Rx = [1.0       0.0      0.0
          0.0  cos(inc) sin(inc)
          0.0 -sin(inc) cos(inc)];
    Rz = [cos(acim) -sin(acim) 0.0
          sin(acim)  cos(acim) 0.0
                0.0        0.0 1.0];
    R = Rx*Rz;
    %Llevo los puntos atraídos al sistema del elipsoide
    ptos = R*(puntos(:,1:3)'-repmat([xe ye ze]',1,nP));
    ptos = ptos';
    %Compruebo el algoritmo a utilizar
    if algoritmo==0
        %Calculo la atracción del elipsoide sobre los puntos
        atrac = GravElipsoideProlatoY([elip(i,4) elip(i,5) elip(i,8)],ptos);
    else
        %Discretizo el elipsoide
        [pelip,aux] = DiscretizaElipsoideProlato(elip(i,4),elip(i,5),pe,tmax);
        n = size(pelip,1);
        %Compruebo los argumentos de salida
        if nargout>2
            %Asigno las diferencias relativas de volumen
            dvr(i) = aux;
            if nargout>3
                %Llevo las coordenadas de los puntos al sistema global
                aux = R'*pelip'+repmat([xe ye ze]',1,n);
                %Asigno las coordenadas
                edisc{i} = aux';
            end
        end
        %Añado la masa de cada punto a la matriz
        m = 4.0/3.0*pi*elip(i,4)*elip(i,5)^2*elip(i,8);
        m = m/n;
        pelip = [pelip m*ones(n,1)];
        %Calculo la atracción del elipsoide sobre los puntos
        atrac = GravNubePuntos(pelip,ptos);
    end
    %Llevo los vectores atracción al sistema global
    atrac = R'*atrac';
    %Actualizo la matriz de atracciones
    g = g+atrac';
    %Compruebo si se pide el segundo argumento de salida
    if nargout>1
        %Genero los puntos sobre la superficie del elipsoide
        [aux,y,z] = ellipsoid(0.0,0.0,0.0,elip(i,5),elip(i,4),elip(i,5));
        aux = [aux(:) y(:) z(:)];
        %Los transformo al sistema global
        aux = R'*aux'+repmat([xe ye ze]',1,size(aux,1));
        %Los asigno al array de salida
        esup{i} = aux';
    end
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
