clear all;
clf;
hold off;

#===================================
#       PLOT FUNCTIONS
#===================================

function plot_object(ax, object)
  V = object.V
  E = object.E

  n_rows = size(E)(1);

  for i=1:1:n_rows
    v1 = E(i, 1);
    v2 = E(i, 2);

    X1 = V(1, v1);
    Y1 = V(2, v1);
    Z1 = V(3, v1);

    X2 = V(1, v2);
    Y2 = V(2, v2);
    Z2 = V(3, v2);

    plot3([Z1 Z2], [X1 X2], [Y1 Y2], 'color', object.color, ax=ax)
    hold on
  endfor
end

function scatter_object(ax, object)
  V = object.V

  [r c] = size(V)

  for i=1:1:c
    X = V(1, i)
    Y = V(2, i)
    Z = V(3, i)

    scatter3(Z, X, Y, 6, "filled", "marker")
    hold on
  endfor
end

#=========================================
#     OBJECT TRANSFORMATION FUNCTIONS
#=========================================

function object = object_scale(old_object, A, B, C)
  M = [
    A 0 0 0 ;
    0 B 0 0 ;
    0 0 C 0 ;
    0 0 0 1
  ]

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  )
end

function object = object_rotate_Y(old_object, angle)
  A = cos(angle)
  B = sin(angle)

  M = [
    +A +0 +B +0 ;
    +0 +1 +0 +0 ;
    -B +0 +A +0 ;
    +0 +0 +0 +1
  ]

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  )
end

function object = object_rotate_X(old_object, angle)
  A = cos(angle)
  B = sin(angle)

  M = [
    +1 +0 +0 +0 ;
    +0 +A -B +0 ;
    +0 +B +A +0 ;
    +0 +0 +0 +1
  ]

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  )
end

function object = object_rotate_Y(old_object, angle)
  A = cos(angle)
  B = sin(angle)

  M = [
    +A +0 +B +0 ;
    +0 +1 +0 +0 ;
    -B +0 +A +0 ;
    +0 +0 +0 +1
  ]

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  )
end

function object = object_rotate_Z(old_object, angle)
  A = cos(angle)
  B = sin(angle)

  M = [
    +A -B +0 +0 ;
    +B +A +0 +0 ;
    -0 +0 +1 +0 ;
    +0 +0 +0 +1
  ]

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  )
end

function object = object_translate(old_object, tx, ty, tz)
  M = [
    1 0 0 tx ;
    0 1 0 ty ;
    0 0 1 tz ;
    0 0 0 1
  ]

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  )
end

#===================================
#     OBJECT CREATION FUNCTIONS
#===================================

function object = create_cube(color)
  object = struct(
  'V', [
    +0.5 -0.5 -0.5 +0.5 +0.5 -0.5 -0.5 +0.5 ;
    +0.5 +0.5 +0.5 +0.5 -0.5 -0.5 -0.5 -0.5 ;
    -0.5 -0.5 +0.5 +0.5 -0.5 -0.5 +0.5 +0.5 ;
    +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0 +1.0
  ],
  'E', [
    1 2 ; 2 3 ; 3 4 ; 4 1 ; 5 6 ; 6 7 ; 7 8 ; 8 5 ; 1 5 ; 2 6 ; 3 7 ; 4 8
  ],
  'color', color
)
end

function object = create_cylinder(
  circles,    # number of circles
  pp_circle,  # points per circle
  color
)
  i = 1:1:pp_circle;
  i = i / pp_circle;
  i = i * 2 * pi;

  circle_X = cos(i);
  circle_Y = -0.5 * ones([1 length(i)]);
  circle_Z = sin(i);

  X = [];
  Y = [];
  Z = [];

  i = 1:1:pp_circle;
  i = i';
  circle_E = [i i+1];
  circle_E(pp_circle, 2) = 1;

  E = [];

  step = 1 / (circles - 1);

  for i=0:(circles-1)
    X = [X circle_X];
    Y = [Y (circle_Y + step * i)];
    Z = [Z circle_Z];
  endfor

  for i=0:(circles-1)
    E = [E ; (circle_E + pp_circle * i)];
  endfor

  i = (0:circles-2)';
  i = [i (i + 1)];
  i = i * 20 + 1;

  for j=0:(pp_circle - 1)
    E = [E ; i + j];
  endfor

  object = struct(
    'V', [
      X ;
      Y ;
      Z ;
      ones([1 pp_circle * circles]) ;
    ],
    'E', E,
    'color', color
  );
end

function object = create_cone(bases, pp_base, color)
  t = 1:bases;
  t = t / bases;

  V_top = [ 0 ; 0.5 ; 0 ; 1];
  V_o= [
    zeros([1 length(t)]) ;
    -t + 0.5 ;
    0.5 * t ;
    ones([1 length(t)])
  ];

  t = 1:pp_base;
  t = t / pp_base;
  t = 2 * pi * t;

  V = []

  for t=t
    s = sin(t);
    c = cos(t);

    M = [
      +c +0 +s +0 ;
      +0 +1 +0 +0 ;
      -s +0 +c +0 ;
      +0 +0 +0 +1 ;
    ];

    V_n = M * V_o;

    V = [V V_n];
  endfor

  V = [V V_top]

  top_index = pp_base * bases + 1

  t = 1:bases:1+bases*(pp_base - 1);
  t = t';

  u = top_index * ones([length(t) 1])

  E = []
  for j=1:bases
    E = [E ; [u t]];
    E = [E ; [t shift(t, length(t)-1)]]
    u = t;

    t = t + 1
  endfor

  object = struct(
    'V', V,
    'E', E,
    'color', color
  );
end

ax = gca()

# cube = create_cube("#00FF00")
# cube = object_scale(cube, 2, 2, 2)
# cube = object_translate(cube, -2.5, 0, 0)
# plot_object(ax, cube)

#cylinder = create_cylinder(20, 20, '#FF0000')
#cylinder = object_scale(cylinder, 1, 2, 1)
# cylinder = object_translate(cylinder, 2.5, 0, 0)
#plot_object(ax, cylinder)

cone = create_cone(20, 20, '#0000FF')
plot_object(ax, cone)

set(ax, "xdir", "reverse")
set(ax, "ydir", "reverse")

set(ax, "xlabel", "Z")
set(ax, "ylabel", "X")
set(ax, "zlabel", "Y")

set(ax, "xlim", [-1 +1])
set(ax, "ylim", [-1 +1])
set(ax, "zlim", [-1 +1])
