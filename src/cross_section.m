% Matboard Bridge - small-scale box girder bridge
% Copyright © 2023 Nabeth Ghazi
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

classdef cross_section
    % cross-section of a box girder bridge made up of vertical/hoizontal rectangles

    properties
        bi
        hi
        yi
        h_total
        y_glue
        % buck 1-3: flexural buckling
        % buck V: shear buckling
        y_buck1
        y_buck2
        y_buck3

        t_buck1
        t_buck2
        t_buck3
        t_buckV

        b_buck1
        b_buck2
        b_buck3
        b_buckV

        a_buckV
    end

    methods

        function obj = cross_section(bi, hi, yi, h_total, y_glue, y_buck1, ...
                y_buck2, y_buck3, t_buck1, t_buck2, t_buck3, t_buckV, ...
                b_buck1, b_buck2, b_buck3, b_buckV, a_buckV)
            obj.bi = bi;
            obj.hi = hi;
            obj.yi = yi;
            obj.h_total = h_total;
            obj.y_glue = y_glue;

            obj.y_buck1 = y_buck1;
            obj.y_buck2 = y_buck2;
            obj.y_buck3 = y_buck3;

            obj.t_buck1 = t_buck1;
            obj.t_buck2 = t_buck2;
            obj.t_buck3 = t_buck3;
            obj.t_buckV = t_buckV;

            obj.b_buck1 = b_buck1;
            obj.b_buck2 = b_buck2;
            obj.b_buck3 = b_buck3;
            obj.b_buckV = b_buckV;

            obj.a_buckV = a_buckV;
        end


        function ybar = ybar(obj)
            Ai = obj.bi .* obj.hi;
            ybar = sum(Ai .* obj.yi) / sum(Ai);
        end


        function ytop = ytop(obj)
            ytop = obj.ybar - obj.h_total;
        end


        function ybot = ybot(obj)
            ybot = obj.ybar - 0;
        end


        function I = I(obj)
            Ai = obj.bi .* obj.hi;
            di = obj.yi - obj.ybar;
            Ii = (obj.bi .* obj.hi.^3) / 12;
            I = sum(Ii + Ai .* di.^2);
        end


        function Q = Q(obj, y)
            ytop = obj.yi + obj.hi / 2;
            ybot = obj.yi - obj.hi / 2;

            if y >= obj.ybar
                b = obj.bi(ytop > y);
                ybot = ybot(ytop > y);
                ytop = ytop(ytop > y);

                ybot(ybot < y) = y;

            else
                b = obj.bi(ybot < y);
                ytop = ytop(ybot < y);
                ybot = ybot(ybot < y);

                ytop(ytop > y) = y;
            end

            h = (ytop - ybot);
            d = abs(ybot + h / 2 - obj.ybar);
            A = b .* h;
            Q = sum (A .* d);

        end


        function bV = bV(obj, y)
            ytop = obj.yi + obj.hi / 2;
            ybot = obj.yi - obj.hi / 2;

            ep = 1e-9;
            bV_top = sum (obj.bi(ybot <= (y + ep) & (y + ep) <= ytop));
            bV_bot = sum (obj.bi(ybot <= (y - ep) & (y - ep) <= ytop));
            bV = min(bV_top, bV_bot);
        end

    end

    methods (Static)

        function S_buck_ult = S_buck_ult(k, E, mu, t, b)
            S_buck_ult = (k * pi^2 * E) / (12 * (1 - mu^2)) * (t / b)^2;
        end


        function T_buck_ult = T_buck_ult(k, E, mu, t, b, a)
            T_buck_ult = (k * pi^2 * E) / ...
                (12 * (1 - mu^2)) * ((t / b)^2 + (t / a)^2);
        end

    end

end