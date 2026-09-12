% Parameters
time_step = 0.0005;  % Time step (s)
total_time = 5.0;  % Total simulation time (s)
g = 9.81;  % Gravitational acceleration (m/s^2)
particle_radius = 0.01;  % Radius of particles (m)
particle_radius_min = 0.01;  % Radius of particles (m)
particle_radius_max = 0.05;  % Radius of particles (m)
particle_mass = 0.1;  % Mass of particles (kg)
num_particles = 100;  % Number of particles
container_width = 1.0;  % Width of the container (m)
container_height = 1.0;  % Height of the container (m)
restitution_coefficient = 0.8;  % Coefficient of restitution for collisions

% Initialize particle positions and velocities
rng(42);  % For reproducibility
positions = rand(num_particles, 2) .* [container_width, container_height/3] + [0, 2*container_height/3];
velocities = zeros(num_particles, 2);
forces = zeros(num_particles, 2);

% Pipe segments
num_segments = 6;
pipe_width = 0.02;
pipe_pos_x = 0.4;
segments = zeros(num_segments, 4);
segments(1, :) = [0.1, 0.7, pipe_pos_x-pipe_width, 0.5];
segments(2, :) = [pipe_pos_x-pipe_width, 0.5, pipe_pos_x-pipe_width, 0.1];
segments(3, :) = [pipe_pos_x+pipe_width, 0.1, pipe_pos_x+pipe_width, 0.5];
segments(4, :) = [pipe_pos_x+pipe_width, 0.5, 0.9, 0.7];
ves_width = 0.2;
ves_height = 0.1;
segments(5, :) = [pipe_pos_x - ves_width/2, ves_height, pipe_pos_x - ves_width/2, 0];
segments(6, :) = [pipe_pos_x + ves_width/2, ves_height, pipe_pos_x + ves_width/2, 0];

% Handle particle-particle collisions
function handle_collisions(positions, velocities, particle_radius, particle_mass, restitution_coefficient)
    for i = 1:num_particles
        for j = i+1:num_particles
            delta_pos = positions(i, :) - positions(j, :);
            distance = norm(delta_pos);
            if distance < 2 * particle_radius
                normal = delta_pos / distance;
                relative_velocity = velocities(i, :) - velocities(j, :);
                relative_speed = dot(relative_velocity, normal);

                if relative_speed < 0
                    impulse = -(1 + restitution_coefficient) * relative_speed / (2 / particle_mass);
                    velocities(i, :) = velocities(i, :) + impulse * normal / particle_mass;
                    velocities(j, :) = velocities(j, :) - impulse * normal / particle_mass;
                end

                overlap = 2 * particle_radius - distance;
                positions(i, :) = positions(i, :) + normal * overlap / 2;
                positions(j, :) = positions(j, :) - normal * overlap / 2;
            end
        end
    end
end

% Handle wall collisions
function handle_wall_collisions(positions, velocities, container_width, container_height, particle_radius, restitution_coefficient)
    for i = 1:num_particles
        if positions(i, 1) - particle_radius < 0 || positions(i, 1) + particle_radius > container_width
            velocities(i, 1) = velocities(i, 1) * -restitution_coefficient;
            if positions(i, 1) - particle_radius < 0
                positions(i, 1) = particle_radius;
            end
            if positions(i, 1) + particle_radius > container_width
                positions(i, 1) = container_width - particle_radius;
            end
        end

        if positions(i, 2) - particle_radius < 0 || positions(i, 2) + particle_radius > container_height
            velocities(i, 2) = velocities(i, 2) * -restitution_coefficient;
            if positions(i, 2) - particle_radius < 0
                positions(i, 2) = particle_radius;
            end
            if positions(i, 2) + particle_radius > container_height
                positions(i, 2) = container_height - particle_radius;
            end
        end
    end
end

% Detect collision with line segments and reflect velocities
function detect_collision_and_reflect(segments, positions, velocities, particle_radius, restitution_coefficient)
    for i = 1:size(segments, 1)
        x1 = segments(i, 1);
        y1 = segments(i, 2);
        x2 = segments(i, 3);
        y2 = segments(i, 4);
        segment_vector = [x2 - x1, y2 - y1];
        for k = 1:num_particles
            xs = positions(k, 1);
            ys = positions(k, 2);
            R = particle_radius;
            line_to_circle = [xs - x1, ys - y1];

            t_closest = dot(line_to_circle, segment_vector) / dot(segment_vector, segment_vector);
            t_closest_clipped = max(0, min(1, t_closest));
            closest_point = [x1, y1] + t_closest_clipped * segment_vector;
            distance = norm(closest_point - [xs, ys]);

            if distance <= R
                normal = (closest_point - [xs, ys]) / distance;
                velocity = velocities(k, :);
                velocities(k, :) = velocity - 2 * dot(velocity, normal) * normal;
                positions(k, :) = closest_point - normal * particle_radius;
            end
        end
    end
end

% Update particle positions
function update_positions(positions, velocities, time_step, g)
    velocities(:, 2) = velocities(:, 2) - g * time_step;
    positions = positions + velocities * time_step;
end

% Animation update function
function update(frame,positions,velocities,time_step,g,particle_radius,particle_mass,restitution_coefficient,num_particles)
    update_positions(positions, velocities, time_step, g);
    handle_collisions(positions, velocities, particle_radius, particle_mass, restitution_coefficient);
    handle_wall_collisions(positions, velocities, container_width, container_height, particle_radius, restitution_coefficient);
    detect_collision_and_reflect(segments, positions, velocities, particle_radius, restitution_coefficient);
    set(scat, 'XData', positions(:, 1), 'YData', positions(:, 2));
end

% Set up the figure and axis
figure;
hold on;
axis([0 container_width 0 container_height]);
scat = scatter(positions(:, 1), positions(:, 2), 2000*particle_radius, 'filled', 'MarkerFaceColor', 'blue');
for i = 1:size(segments, 1)
    plot([segments(i, 1), segments(i, 3)], [segments(i, 2), segments(i, 4)], 'g-', 'LineWidth', 2);
end

% Create animation
while true
    update(0,positions,velocities,time_step,g,particle_radius,particle_mass,restitution_coefficient,num_particles);
    pause(time_step);
end
