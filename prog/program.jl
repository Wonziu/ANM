using Plots
using printf

# WSZYSTKIE FUNKCJE, KTÓRYCH UŻYŁEM DO WYKONANIA OBLICZEŃ

# Funkcja która zwraca prezenetacje działania metody bisekcji
function getbisection()

    f(x) = x^2 - 3

    default(tickfont = (12, :orange), framestyle = :zerolines)
    plot(f, 0, 6, linewidth = 2, label = "f(x) = x^2 - 3", legend = :topleft)
    annotate!((1, -2, "a"), (5, -2, "b"), (sqrt(3), -2, "alfa"), (3, -2, "c"))
    scatter!([(1, 0), (5, 0), (sqrt(3), 0), (3, 0)], c = :red, label = false)

    savefig("bisectionPlot.pdf")
end

# Funkcja która zwraca prezenetacje działania metody Newtona
function getNewton()
    x0 = 5

    f(x) = x^2 - 3
    g(x) = 2x

    h(x) = g(x0) * x + f(x0) - g(x0) * x0

    default(tickfont = (12, :orange), framestyle = :zerolines)
    plot([f, h], 0, 7, linewidth = 2, label = ["f(x) = x^2 - 3" "styczna w x0"], legend = :topleft)
    annotate!((5.6, 22, "f(x0)"), (2.8, -5, "x1"), (5, -5, "x0"), (1.7, -5, "alfa"))
    scatter!([(5, 22), (5, 0), (2.8, 0), (sqrt(3), 0)], c = :red, label = false)

    savefig("newtonPlot.pdf")
end

# Funkcja która wartość funkcji w punkcie x, oraz jej pochodnej w punkcie x
function Horner(polynomial, x)
    n = length(polynomial) - 1
    
    f  = polynomial[n]
    ∂f = 0
    
    for i = n - 1:-1:0
        ∂f = f + x * ∂f
        f  = polynomial[i] + x * f
    end
    
    return f, ∂f
end

# Funkcja dzieli wielomian 3 stopnia, przez wielomian (x - x0)
function divide(f, x0) 
    b = zeros(3)
    b[3] = f[3]
    
    for i = 2:-1:1
        b[i] = f[i] + x0 * b[i + 1]
    end
 
    return Polynomial(b)
end

#Funkcja znajduje pieriawstki równania kwadratowego
function solveQuadratic(f)
    a = f[2]
    b = f[1]
    c = f[0]
    
    x1 = 0
    x2 = 0
    
    Δ = b * b - 4 * a * c
    
    if (Δ > 0)
        
        if (b > 0)
            x1 = (-b - sqrt(Δ)) / (2 * a)
        else
            x1 = (-b + sqrt(Δ)) / (2 * a)
        end
        x2 = c / (a * x1)
            
    elseif (Δ == 0)        
        x1 = -b / (2 * a)        
        x2 = -b / (2 * a)
    else
        
        if (b > 0)
            x1 = (-b - sqrt(-Δ)*im) / (2 * a)
        else
            x1 = (-b + sqrt(-Δ)*im) / (2 * a)
        end
        x2 = c / (a * x1)
    end
    
    return [x1, x2]
end

# Zaimplementowana metoda Newtona
function Newton(f, x0, ε, M, print = false)
    v = Horner(f, x0)
    
    if (abs(v[1]) < ε)
        return x0
    end
    
    for i = 1:M
        x1 = x0 - v[1] / v[2]
        v = Horner(f, x1)
        
        rε = abs((x1 - x0) / x1)
        
        if print
            @printf("Iteracja: %d, x_i = %f, błąd bezwględny = %f\n", i, x1, rε)
        end
            
        if (abs(v[1]) < ε)
            x0, x2 = solveQuadratic(divide(f, x1))
            
            return [x0, x1, x2, f(x1)]
        end
        
        x0 = x1
    end
    
    x1, x2 = solveQuadratic(divide(f, x0))
            
    return [x0, x1, x2, f(x0)]
end

# Zaimplementowana Metoda Bisekcji
function Bisection(f, x0, x1, ε, M, print = false)
    a = Horner(f, x0)[1]
    b = Horner(f, x1)[1]
    
    if (sign(a) == sign(b))
        x1, x2 = solveQuadratic(divide(f, x0))
                
        return [x0, x1, x2]
    end
    
    for i = 1:M
        x2 = ((x0 + x1) / 2)
        
        v = Horner(f, x2)[1]
        
        if print
            @printf("Iteracja: %d, x0 = %f, f(x0) = %f\n", i, x2, v)
        end
        
        if (abs(v) < ε)
            x0, x1 = solveQuadratic(divide(f, x2))
            
            return [x0, x1, x2, f(x0)]
        end
        
        if (sign(a) != sign(v))
            x1 = x2
        else
            x0 = x2
        end
    end
        
    x1, x2 = solveQuadratic(divide(f, x0))
                
    return [x0, x1, x2, f(x0)]
end

# Zaimplementowana metoda Bairstowa
function Bairstow(f, u, v, M)
    n = length(f)
    
    b = zeros(n)
    c = zeros(n)
    
    b[n] = f[n - 1]
    c[n - 1] = f[n - 1]
    
    for j = 1:M
        b[n-1] = f[n - 2] + u * b[n]
        for k = n - 2:-1:1
            b[k] = f[k - 1] + u * b[k + 1] + v * b[k + 2]
            c[k] = b[k + 1] + u * c[k + 1] + v * c[k + 2]
        end
        
        J = c[1] * c[3] - c[2] * c[2]
        u = u + (c[2] * b[2] - c[3] * b[1]) / J
        v = v + (c[2] * b[1] - c[1] * b[2]) / J
    end
    
    quadratic = Polynomial([-v, -u, 1])
    x0, x1 = solveQuadratic(quadratic)
    
    if (x0 != 0 && x1 != 0)
        x2 = -f[0] / (f[3] * x0 * x1)
    else
        x2 = -f[2] / f[3] - x0 - x1
    end
    
    return [x0, x1, x2, f(x0)]
    
end

# Zaimplementowane wzory Cardano
function Cardano(f)
    a = f[3]
    b = f[2]
    c = f[1]
    d = f[0]
    
    Q = (3 * a * c - b * b) / (9 * a * a)
    R = (9 * a * b * c - 27 * a * a * d - 2 * b * b * b) / (54 * a * a * a)
    
    discriminant = Q * Q * Q + R * R
    S = 0
    T = 0
    
    if (discriminant >= 0)
        S = cbrt(R + sqrt(discriminant))
        T = cbrt(R - sqrt(discriminant))
    else
        S = (R + sqrt(discriminant + 0im))^(1/3)
        T = (R - sqrt(discriminant + 0im))^(1/3)
    end
    
    x0 = S + T - b / (3 * a)
    x1 = -(S + T) / 2 - b / (3 * a) + (im * sqrt(3)) / 2 * (S - T)
    x2 = -(S + T) / 2 - b / (3 * a) - (im * sqrt(3)) / 2 * (S - T)
    
    return [x0, x1, x2]
end

# Funkcja która zwraca wykres przedstawiający liczbę cyfr znaczących dla każdego pierwiastka dla metody Bairstowa
function TestBairstow(f, v, u, MaxM, r1, r2, r3)
    
    function significantDigits(x)
        if (x == 0)
            return 20
        end
        
        return min(20, max(0, round(-log10(abs(x)))))
    end
    
    x0 = zeros(MaxM)
    x1 = zeros(MaxM)
    x2 = zeros(MaxM)
    
    for i = 1:MaxM
        b = Bairstow(f, u, v, i)
    
        x0[i] = significantDigits(abs(r1 - b[1]))
        x1[i] = significantDigits(abs(r2 - b[2]))
        x2[i] = significantDigits(abs(r3 - b[3]))        
    end    
    
    default(framestyle = :zerolines)
    plot([x0, x1, x2], label = ["x0" "x1" "x2"], 
        linewidth = 3, xlabel = "Iteracja", ylabel="Liczba cyfr znaczących", legend = :topleft)    
end

# Funkcja, która znajduje wstępne przybliżenie pierwiastka dla metody Newtona
function getInitialGuess(f, df)
    a = f[3]
    
    if (typeof(solveQuadratic(df)[1])<:Real)
    
        d0 = sort(solveQuadratic(df))

        f1 = Horner(f, d0[1] - 0.5)[1]
        f2 = Horner(f, d0[2] + 0.5)[1]

        if (sign(a) == sign(f1))
            return d0[1] - 0.5
        elseif (sign(a) != sign(f2))
            return d0[2] + 0.5
        else
            return d0[2] + 0.1
        end
    end
    
    return -0.350001
end

# Funkcja która zwraca wykres przedstawiający liczbę cyfr znaczących dla każdego pierwiastka dla metody Newtona
function TestNewton(f, x_0, ε, MaxM, r1, r2, r3)
    
    function significantDigits(x)
        if (x == 0)
            return 20
        end
        
        return min(20, max(0, round(-log10(abs(x)))))
    end
    
    x0 = zeros(MaxM)
    x1 = zeros(MaxM)
    x2 = zeros(MaxM)
    
    for i = 1:MaxM
        b = Newton(f, x_0, ε, i)
    
        x0[i] = significantDigits(abs(r1 - b[1]))
        x1[i] = significantDigits(abs(r2 - b[2]))
        x2[i] = significantDigits(abs(r3 - b[3]))        
    end    
    
    default(framestyle = :zerolines)
    plot([x0, x1, x2], label = ["x0" "x1" "x2"], 
        linewidth = 3, xlabel = "Iteracja", ylabel="Liczba cyfr znaczących", legend = :topleft)    
end

# Funkcja która zwraca wykres przedstawiający liczbę cyfr znaczących dla każdego pierwiastka dla metody bisekcji
function TestBisection(f, x_0, x_1, ε, MaxM, r1, r2, r3)
    
    function significantDigits(x)
        if (x == 0)
            return 20
        end
        
        return min(20, max(0, round(-log10(abs(x)))))
    end
    
    x0 = zeros(MaxM)
    x1 = zeros(MaxM)
    x2 = zeros(MaxM)
    
    for i = 1:MaxM
        b = Bisection(f, x_0, x_1, ε, i)

        x0[i] = significantDigits(abs(r1 - b[1]))
        x1[i] = significantDigits(abs(r2 - b[2]))
        x2[i] = significantDigits(abs(r3 - b[3])) 
    end    
    
    default(framestyle = :zerolines)
    plot([x0, x1, x2], label = ["x0" "x1" "x2"], 
        linewidth = 3, xlabel = "Iteracja", ylabel="Liczba cyfr znaczących", legend = :topleft)    
end

# Funkcja która zwraca wartości bezwzględne pierwiastków korzystając ze wzorów Cardano
function TestCardano(f, r1, r2, r3)
    
    x = Cardano(f)
    
    println("Błąd względny x0: ", abs(x[1] - r1))
    println("Błąd względny x1: ", abs(x[2] - r2))
    println("Błąd względny x2: ", abs(x[3] - r3)) 
end