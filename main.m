clear all;
clf;
hold off;
system('clear');

#===================================
#       PLOT FUNCTIONS
#===================================

function plot_object(object)
  V = object.V;
  E = object.E;

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

    plot3([Z1 Z2], [X1 X2], [Y1 Y2], 'color', object.color);
    hold on;
  endfor
end

function scatter_object(object)
  V = object.V;

  [r c] = size(V);

  for i=1:1:c
    X = V(1, i);
    Y = V(2, i);
    Z = V(3, i);

    scatter3(Z, X, Y, 6, [255, 0, 0], "filled");
    hold on;
  endfor
end

function plot_2d_object(object, M)
  V = object.V;
  E = object.E;

  V = M * V;

  n_edges = size(E)(1);
  c = size(V)(2);
  t = 1:c;

  for x=t
    V(1, x) = V(1, x) / V(4, x);
    V(2, x) = V(2, x) / V(4, x);
    V(3, x) = V(3, x) / V(4, x);
  endfor

  for x=1:n_edges
    V1 = E(x, 1);
    V2 = E(x, 2);

    X1 = V(1, V1);
    Y1 = V(2, V1);

    X2 = V(1, V2);
    Y2 = V(2, V2);

    plot([X1 X2], [Y1 Y2], 'color', object.color);
    hold on;
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
  ];

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  );
end

function object = object_rotate_X(old_object, phi)
  A = cos(phi);
  B = sin(phi);

  M = [
    +1 +0 +0 +0 ;
    +0 +A -B +0 ;
    +0 +B +A +0 ;
    +0 +0 +0 +1
  ];

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  );
end

function object = object_rotate_Y(old_object, phi)
  A = cos(phi);
  B = sin(phi);

  M = [
    +A +0 +B +0 ;
    +0 +1 +0 +0 ;
    -B +0 +A +0 ;
    +0 +0 +0 +1
  ];

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  );
end

function object = object_rotate_Z(old_object, phi)
  A = cos(phi);
  B = sin(phi);

  M = [
    +A -B +0 +0 ;
    +B +A +0 +0 ;
    -0 +0 +1 +0 ;
    +0 +0 +0 +1
  ];

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  );
end

function object = object_translate(old_object, tx, ty, tz)
  M = [
    1 0 0 tx ;
    0 1 0 ty ;
    0 0 1 tz ;
    0 0 0 1
  ];

  object = struct(
    'V', M * old_object.V,
    'E', old_object.E,
    'color', old_object.color
  );
end

function M = create_camera_system_transformation_matrix(at, eye)
  n = (at - eye)(1:3, 1);
  n = n / sqrt(n' * n);

  aux = [ 0 ; 1 ; 0 ];
  proj = (aux' * n) / (n' * n) * n;

  v = aux - proj;
  v = v / sqrt(v' * v);

  u = cross(v, n);

  R = [
    u' 0 ;
    v' 0 ;
    n' 0 ;
    0 0 0 1;
  ]

  T = [
    1 0 0 -eye(1, 1) ;
    0 1 0 -eye(2, 1) ;
    0 0 1 -eye(3, 1) ;
    0 0 0 1 ;
  ]

  M = R * T;
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
  );
end

function object = create_cylinder(
  circles,    # number of circles
  pp_circle,  # points per circle
  color
)
  t = 0:circles-1;
  t = t / (circles - 1);

  V_o = [
    zeros([1 circles]) ;
    0.5 - t ;
    zeros([1 circles]) + 1;
    ones([1 circles])
  ];

  V = [];

  t = 0:pp_circle-1;
  t = t / pp_circle;
  t = t * 2 * pi;

  for x=t
    A = cos(x);
    B = sin(x);

    M = [
      +A +0 +B +0 ;
      +0 +1 +0 +0 ;
      -B +0 +A +0 ;
      +0 +0 +0 +1 ;
    ];

    V_n = M * V_o;
    V = [V V_n];
  endfor

  n_points = pp_circle * circles;
  t = 1:n_points;

  t = reshape(t, circles, pp_circle);
  u = shift(t, -1, 1);

  t = t(1:circles-1, :);
  u = u(1:circles-1, :);

  E = [reshape(t, n_points - pp_circle, 1) reshape(u, n_points - pp_circle, 1)];

  t = 1:n_points;

  t = reshape(t, circles, pp_circle)';
  u = shift(t, -1, 1);

  E = [E ; [reshape(t, n_points, 1) reshape(u, n_points, 1)]];


  object = struct(
    'V', V,
    'E', E,
    'color', color
  );
end

function object = create_cone(bases, pp_base, color)
  t = 1:bases;
  t = t / bases;

  V_top = [ 0 ; 0.5 ; 0 ; 1 ];

  V_o= [
    zeros([1 length(t)]) ;
    -t + 0.5 ;
    t ;
    ones([1 length(t)])
  ];

  t = 1:pp_base;
  t = t / pp_base;
  t = 2 * pi * t;

  V = [];

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

  V = [V V_top];

  top_index = pp_base*bases + 1;

  t = 1:top_index-1;
  t = reshape(t, bases, pp_base);
  t = [zeros([1 pp_base]) + top_index ; t];

  s = shift(t, -1, 1);

  t = t(1:bases, :);
  s = s(1:bases, :);

  t = reshape(t, [top_index-1 1]);
  s = reshape(s, [top_index-1 1]);

  E = [ t s ];

  t = 1:top_index-1;
  t = reshape(t, bases, pp_base)';

  s = shift(t, -1, 1);

  t = reshape(t, [top_index-1 1]);
  s = reshape(s, [top_index-1 1]);

  E = [E ; [ t s ]];

  object = struct(
    'V', V,
    'E', E,
    'color', color
  );
end

# Creates a sphere of radius equal to r
# @param r          Radius of the sphere.
# @param n_arcs  Amount o arcs that form the sphere (>= 2)
# @param pp_circle  Amount of points that will form each arc (>= 3)
function object = create_sphere(
  n_arcs,
  pp_arc,
  color
)
  t = (1:pp_arc-2) / (pp_arc - 1);
  t = (t * pi) - (pi / 2);

  V_arc = [ cos(t) ; sin(t) ; zeros([1 pp_arc-2]) ; ones([1 pp_arc-2]) ];

  t = 0:n_arcs-1;
  t = t / n_arcs;
  t = 2 * pi * t;

  V = [];
  for x=t
    A = cos(x);
    B = sin(x);
    M = [
      +A +0 +B +0 ;
      +0 +1 +0 +0 ;
      -B +0 +A +0 ;
      +0 +0 +0 +1 ;
    ];

    V_new_arc = M * V_arc;
    V = [V V_new_arc];
  endfor

  V_top = [ 0 ; 1 ; 0 ; 1 ];
  V_bottom = [ 0 ; -1 ; 0 ; 1 ];

  V = [V V_top V_bottom ];

  n_points = (pp_arc-2)*n_arcs + 2;

  t = 1:n_points-2;
  t = reshape(t, pp_arc-2, n_arcs);

  t = cat(1, n_points + zeros([1 n_arcs]), t, n_points - 1 + zeros([1 n_arcs]));

  u = shift(t, -1, 1);

  t = t(1:pp_arc-1, 1:n_arcs);
  u = u(1:pp_arc-1, 1:n_arcs);

  t = reshape(t, (pp_arc-1)*n_arcs, 1);
  u = reshape(u, (pp_arc-1)*n_arcs, 1);

  E = [t u];

  t = 1:n_points-2;
  t = reshape(t, pp_arc-2, n_arcs)';

  u = shift(t, -1, 1);

  t = reshape(t, (pp_arc-2)*n_arcs, 1);
  u = reshape(u, (pp_arc-2)*n_arcs, 1);

  E = [E ; [t u]];

  object = struct(
    'V', V,
    'E', E,
    'color', color
  );
end

function object = create_toroid(r, R, circles, pp_circle, color)
  t = (0:pp_circle-1) / pp_circle;
  t = 2 * pi * t;

  V_circle = [
    (cos(t) * r) + R;
    (sin(t) * r);
    zeros([1 pp_circle]) ;
    ones([1 pp_circle]) ;
  ];

  t = (0:circles-1) / circles;
  t = 2 * pi * t;

  V = [];
  for x=t
    A = cos(x);
    B = sin(x);

    M = [
      +A +0 +B +0 ;
      +0 +1 +0 +0 ;
      -B +0 +A +0 ;
      +0 +0 +0 +1 ;
    ];

    V_new_circle = M * V_circle;
    V = [ V V_new_circle ];
  endfor

  t = 1:pp_circle*circles;
  t = reshape(t, pp_circle, circles);
  u = shift(t, -1, 1);

  t = reshape(t, pp_circle*circles, 1);
  u = reshape(u, pp_circle*circles, 1);

  E = [t u];

  t = 1:pp_circle*circles;
  t = reshape(t, pp_circle, circles)';
  u = shift(t, -1, 1);

  t = reshape(t, pp_circle*circles, 1);
  u = reshape(u, pp_circle*circles, 1);

  E = [E ; [t u]];

  object = struct(
    'V', V,
    'E', E,
    'color', color
  );
end

function object = create_prysm(k, sides, surfaces, color)
  R = 1 / (2 * sin(pi / sides));
  r = k / (2 * sin(pi / sides));

  t = (0:surfaces-1) / (surfaces-1);

  V_line = [
    zeros([1 surfaces]) ;
    t - 0.5 ;
    (r - R)*t + R;
    ones([1 surfaces])
  ];

  t = (0:sides-1) / sides;
  t = 2 * pi * t;

  V = []
  for x=t
    A = cos(x);
    B = sin(x);

    M = [
      +A +0 +B +0 ;
      +0 +1 +0 +0 ;
      -B +0 +A +0 ;
      +0 +0 +0 +1 ;
    ];

    V_new_line = M * V_line;
    V = cat(2, V, V_new_line);
  endfor

  n_points = sides * surfaces;

  t = 1:n_points;
  t = reshape(t, surfaces, sides);
  u = shift(t, -1, 1);

  t = t(1:surfaces-1, :);
  u = u(1:surfaces-1, :);

  t = reshape(t, n_points - sides, 1);
  u = reshape(u, n_points - sides, 1);

  E = [t u];

  t = 1:n_points;
  t = reshape(t, surfaces, sides)';
  u = shift(t, -1, 1);

  t = reshape(t, n_points, 1);
  u = reshape(u, n_points, 1);

  E = cat(1, E, [t u]);

  object = struct(
    'V', V,
    'E', E,
    'color', color
  );
end

function configure_axis(xlimits, ylimits, zlimits)
  ax = gca()

  set(ax, "xdir", "reverse");
  set(ax, "ydir", "reverse");

  set(ax, 'cameraviewanglemode', 'manual')

  xlabel('Z');
  ylabel('X');
  zlabel('Y');

  axis([ zlimits xlimits ylimits ]);
end

# Creation of the world objects
cube = create_cube('#00FF00');
cube = object_scale(cube, 2, 2, 2);

sphere = create_sphere(20, 20, '#9933FF');
sphere = object_scale(sphere, 1, 1, 1);

cylinder = create_cylinder(4, 20, '#00CCCC');
cylinder = object_scale(cylinder, 1, 3, 1);

cone = create_cone(4, 20, '#0000FF');
cone = object_scale(cone, 1, 3, 1);

frustum = create_prysm(0.6, 4, 4, '#FF0000');
frustum = object_scale(frustum, 2, 3, 2);

toroid = create_toroid(0.5, 1, 20, 10, '#FF8000');
toroid = object_scale(toroid, 2, 2, 2);

# plot_object(cube);
# configure_axis([-1 1], [-1 1], [-1 1]);
# input('Press any key to continue...');
# hold off;

# plot_object(sphere);
# configure_axis([-1 1], [-1 1], [-1 1]);
# input('Press any key to continue...');
# hold off;

# plot_object(cylinder);
# configure_axis([-1.5 1.5], [-1.5 1.5], [-1.5 1.5]);
# input('Press any key to continue...');
# hold off;

# plot_object(cone);
# configure_axis([-1.5 1.5], [-1.5 1.5], [-1.5 1.5]);
# input('Press any key to continue...');
# hold off;

# plot_object(frustum);
# configure_axis([-1.5 1.5], [-1.5 1.5], [-1.5 1.5]);
# input('Press any key to continue...');
# hold off;

# plot_object(toroid);
# configure_axis([-3 3], [-3 3], [-3 3]);
# input('Press any key to continue...');
# hold off;

# Placing the objects in their positions in the scene
cube = object_translate(cube, -2, 1, -3);
cylinder = object_translate(cylinder, -4, +1.5, -7.24);
sphere = object_translate(sphere, -6, 1, -3);
cone = object_translate(cone, 2, +1.5, -3);

frustum = object_rotate_Y(frustum, pi / 4)
frustum = object_translate(frustum, 6, +1.5, -3);

toroid = object_rotate_X(toroid, pi / 6);
toroid = object_translate(toroid, 4, 6, -7.24);

eye = [ -4 ; 5 ; 2.5 ; 1 ];
at = cat(2, cube.V, sphere.V, cylinder.V, cone.V, frustum.V, toroid.V);
at = sum(at, 2) / size(at)(2);

# plot_object(cube);
# plot_object(sphere);
# plot_object(cylinder);
# plot_object(frustum);
# plot_object(cone);
# plot_object(toroid);
# configure_axis([-10 10], [-10 10], [-10 10]);

# input('Press any key to continue...');
# hold off;

# Composing the camera coordinate system
M = create_camera_system_transformation_matrix(at, eye);


origin = M * [ 0 ; 0 ; 0 ; 1];
cube.V = M * cube.V;
cylinder.V = M * cylinder.V;
sphere.V = M * sphere.V;
cone.V = M * cone.V;
frustum.V = M * frustum.V;
toroid.V = M * toroid.V;

plot_object(cube);
plot_object(sphere);
plot_object(cylinder);
plot_object(frustum);
plot_object(cone);
plot_object(toroid);
scatter3(origin(3), origin(1), origin(2), "filled")
configure_axis([-10 10], [-10 10], [0 20]);

input('Press any key to continue');
hold off;

# Projecting the objects in the viewing volume in a bidimensional plane
fovy = 2 * pi / 3;
aspect = 1;
z_far = 20;
top = 1;

z_near = top / tan(fovy / 2)
right = top * aspect

M = [
  z_near / right, 0, 0, 0 ;
  0, -z_near / top, 0, 0 ;
  0, 0, -(z_far + z_near) / (z_far - z_near), -(2 * z_far * z_near) / (z_far - z_near);
  0, 0, -1, 0 ;
]

plot_2d_object(cube, M);
plot_2d_object(cylinder, M);
plot_2d_object(sphere, M);
plot_2d_object(cone, M);
plot_2d_object(frustum, M);
plot_2d_object(toroid, M);
xlim([-right right])
ylim([-1 1]);
