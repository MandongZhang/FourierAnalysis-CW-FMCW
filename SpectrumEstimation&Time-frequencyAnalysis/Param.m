classdef    Param < handle
        properties
            Fs
            Ns
        end
        
        properties(Dependent)
            T
            n
            t
            f
        end
        

        methods
            function obj = Param(Fs,Ns)
                obj.Fs = Fs;
                obj.Ns = Ns;
            end

            function T = get.T(obj)
                T = obj.Ns/obj.Fs;
               
            end

            function n = get.n(obj)
                n = 1:obj.Ns;
            end

            function t = get.t(obj)
                t = obj.T*(obj.n-1)/obj.Ns;
            end
            function f = get.f(obj)
                f = obj.Fs*(obj.n-1)/obj.Ns;
            end
       
        end
end
            