# Test OriHinge
wipe
model BasicBuilder -ndm 3 -ndf 6

# --- Nodos ---
set L 1.0
set theta_0 210
set phi [expr {2.0 * 3.141592653589793 - $theta_0 * 3.141592653589793 / 180.0}]
set x_coord [expr {sin(3.141592653589793/3.0) * $L}]
set z_coord [expr {sin($phi) * $x_coord}]
set x_coord2 [expr {cos($phi) * $x_coord}]

# Ahora los nodos manualmente
node 1 0.0 [expr {-0.5*$L}] 0.0
node 2 0.0 [expr {0.5*$L}] 0.0
node 3 $x_coord 0.0 0.0
node 4 $x_coord2 0.0 $z_coord

# --- Apoyo ---
fix 1 1 1 1 1 1 1
fix 2 1 1 1 1 1 1
fix 3 1 1 1 1 1 1
fix 4 0 0 0 1 1 1

uniaxialMaterial Elastic 1 1000.0

# --- Elementos ---
element corotTruss 1 1 2 1.0 1
element corotTruss 2 2 3 1.0 1
element corotTruss 3 3 1 1.0 1
element corotTruss 4 1 4 1.0 1
element corotTruss 5 2 4 1.0 1
element OriHinge 6 3 1 2 4 0.3

# --- Carga nodal ---
timeSeries Linear 1
pattern Plain 1 1 {
    load 4 0.0 0.0 1 0.0 0.0 0.0
}

# --- Análisis estático ---

system Umfpack
constraints Plain
numberer RCM
test NormDispIncr 1.0e-8 100
algorithm Newton
integrator ArcLength 0.1 1.0
# integrator LoadControl 0.1
analysis Static

# Realizar análisis M veces
puts "Step,Load_Factor,Angle(degrees),Disp_Node4_Z"
set M 10
for {set i 1} {$i <= $M} {incr i} {
    if {[analyze 1] != 0} {
        puts "Análisis fallido en paso $i"
        exit
    }
    set lambda [getLoadFactor 1]
    set theta [eleResponse 6 theta]
    puts "$i,$lambda,[expr {$theta * 180.0 / 3.141592653589793}],[nodeDisp 4 3]" 
}


