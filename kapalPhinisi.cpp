#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iterator>

#define GLM_ENABLE_EXPERIMENTAL

#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// DEKLARASI FUNGSI (FORWARD DECLARATION)
void generateSailMesh(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, int divisions, std::vector<glm::vec3>& target_vertices, std::vector<unsigned int>& target_indices);
float getDeckYAtZ(float z);
float getDeckWidthAtZ(float z);

// VARIABEL GLOBAL & KONFIGURASI

std::vector<glm::vec3> vertices;
std::vector<glm::vec3> normals;
std::vector<unsigned int> indices;

std::vector<glm::vec3> sail_vertices;
std::vector<glm::vec3> sail_normals;
std::vector<unsigned int> sail_indices;

float camera_distance = 80.0f;
float camera_angle_x = 30.0f;
float camera_angle_y = -60.0f;
glm::vec3 camera_target = {0.0f, 15.0f, 0.0f};
bool is_mouse_dragging = false;
int last_mouse_x = 0;
int last_mouse_y = 0;

GLfloat wood_ambient[] = {0.25f, 0.18f, 0.12f, 1.0f};
GLfloat wood_diffuse[] = {0.4f, 0.28f, 0.2f, 1.0f};
GLfloat wood_specular[] = {0.05f, 0.05f, 0.05f, 1.0f};
GLfloat wood_shininess[] = {10.0f};

GLfloat sail_ambient[] = {0.4f, 0.35f, 0.28f, 1.0f};
GLfloat sail_diffuse[] = {0.8f, 0.7f, 0.55f, 1.0f};
GLfloat sail_specular[] = {0.1f, 0.1f, 0.1f, 1.0f};
GLfloat sail_shininess[] = {20.0f};

// STRUKTUR & BLUEPRINT DASAR

struct MastGeometry {
    glm::vec3 base;
    glm::vec3 top;
    glm::vec3 boom_base;
    glm::vec3 boom_tip;
    glm::vec3 gaff_base;
    glm::vec3 gaff_tip;
    glm::vec3 jib_boom_base; // <-- TAMBAHKAN INI
    glm::vec3 jib_boom_tip;  // <-- TAMBAHKAN INI
};

struct HullCrossSection
{
    glm::vec2 P0, P1, P2, P3;
    glm::vec2 getPoint(float t) const {
        float u = 1.0f - t;
        float tt = t * t;
        float uu = u * u;
        float uuu = uu * u;
        float ttt = tt * t;
        return uuu * P0 + 3.0f * uu * t * P1 + 3.0f * u * tt * P2 + ttt * P3;
    }
};

HullCrossSection sternSection =   { {0.0f, 0.0f}, {4.0f, -2.0f}, {5.0f, 9.5f},  {4.0f, 9.0f} };
HullCrossSection aftMidSection =  { {0.0f, 0.0f}, {5.5f, -2.5f}, {7.0f, 8.5f},  {6.5f, 8.8f} };
HullCrossSection midSection =     { {0.0f, 0.0f}, {6.0f, -2.0f}, {7.5f, 8.0f},  {7.0f, 9.0f} };
HullCrossSection foreMidSection = { {0.0f, 0.0f}, {4.0f, -1.5f}, {6.0f, 9.0f},  {5.5f, 9.5f} };
HullCrossSection bowSection =     { {0.0f, 0.0f}, {1.5f, 1.0f},  {4.0f, 9.8f},  {3.0f, 10.0f} };
std::vector<HullCrossSection> blueprint_sections = {sternSection, aftMidSection, midSection, foreMidSection, bowSection};
std::vector<float> blueprint_z_positions = {-22.0f, -15.0f, 0.0f, 15.0f, 25.0f};


// FUNGSI-FUNGSI GENERATOR GEOMETRI

void recalculateNormals(std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref, std::vector<glm::vec3>& norm_ref) {
    norm_ref.assign(vert_ref.size(), glm::vec3(0.0f, 0.0f, 0.0f));
    for (size_t i = 0; i < ind_ref.size(); i += 3) {
        if (ind_ref[i] >= vert_ref.size() || ind_ref[i+1] >= vert_ref.size() || ind_ref[i+2] >= vert_ref.size()) continue;
        glm::vec3 p1 = vert_ref[ind_ref[i]];
        glm::vec3 p2 = vert_ref[ind_ref[i+1]];
        glm::vec3 p3 = vert_ref[ind_ref[i+2]];
        glm::vec3 v1 = p2 - p1;
        glm::vec3 v2 = p3 - p1;
        glm::vec3 normal = glm::cross(v1, v2);
        norm_ref[ind_ref[i]] += normal;
        norm_ref[ind_ref[i+1]] += normal;
        norm_ref[ind_ref[i+2]] += normal;
    }
    for (size_t i = 0; i < norm_ref.size(); ++i) {
        if (glm::length(norm_ref[i]) > 0.0001f) {
            norm_ref[i] = glm::normalize(norm_ref[i]);
        }
    }
}

void generateBox(glm::vec3 pos, glm::vec3 size, std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref) {
    int base_index = vert_ref.size();
    vert_ref.push_back(pos);
    vert_ref.push_back(pos + glm::vec3(size.x, 0, 0));
    vert_ref.push_back(pos + glm::vec3(size.x, size.y, 0));
    vert_ref.push_back(pos + glm::vec3(0, size.y, 0));
    vert_ref.push_back(pos + glm::vec3(0, 0, size.z));
    vert_ref.push_back(pos + glm::vec3(size.x, 0, size.z));
    vert_ref.push_back(pos + glm::vec3(size.x, size.y, size.z));
    vert_ref.push_back(pos + glm::vec3(0, size.y, size.z));
    
    int indices_data[] = {0,1,2, 0,2,3, 1,5,6, 1,6,2, 4,0,3, 4,3,7, 5,4,7, 5,7,6, 3,2,6, 3,6,7, 4,5,1, 4,1,0};
    for(int i=0; i<36; ++i) ind_ref.push_back(base_index + indices_data[i]);
}

void generateCylinder(glm::vec3 pos, float height, float radius, int sides,
                      std::vector<glm::vec3>& vert_ref,
                      std::vector<unsigned int>& ind_ref,
                      const glm::vec3& orientation = {0.0f, 1.0f, 0.0f})
{
    glm::vec3 w = glm::normalize(orientation);
    glm::vec3 start_pos = pos;
    glm::vec3 end_pos = pos + w * height;

    int base_index = vert_ref.size();
    glm::vec3 u, v;
    if (std::abs(w.x) > 0.99f) { 
        u = glm::normalize(glm::cross(w, {0.0f, 1.0f, 0.0f}));
    } else { 
        u = glm::normalize(glm::cross(w, {1.0f, 0.0f, 0.0f})); 
    }
    v = glm::cross(w, u);

    for (int i = 0; i <= sides; ++i) {
        float angle = 2.0f * M_PI * i / sides;
        glm::vec3 offset = radius * (cosf(angle) * u + sinf(angle) * v);
        vert_ref.push_back(start_pos + offset);
    }
    for (int i = 0; i <= sides; ++i) {
        float angle = 2.0f * M_PI * i / sides;
        glm::vec3 offset = radius * (cosf(angle) * u + sinf(angle) * v);
        vert_ref.push_back(end_pos + offset);
    }

    int ring1_start = base_index;
    int ring2_start = base_index + (sides + 1);
    for (int i = 0; i < sides; ++i) {
        ind_ref.push_back(ring1_start + i); ind_ref.push_back(ring2_start + i); ind_ref.push_back(ring1_start + i + 1);
        ind_ref.push_back(ring2_start + i); ind_ref.push_back(ring2_start + i + 1); ind_ref.push_back(ring1_start + i + 1);
    }

    int start_center_idx = vert_ref.size(); vert_ref.push_back(start_pos);
    int end_center_idx = vert_ref.size(); vert_ref.push_back(end_pos);
    for (int i = 0; i < sides; i++) {
        ind_ref.push_back(start_center_idx); ind_ref.push_back(ring1_start + i); ind_ref.push_back(ring1_start + i + 1);
        ind_ref.push_back(end_center_idx); ind_ref.push_back(ring2_start + i + 1); ind_ref.push_back(ring2_start + i);
    }
}

void generateRailSegment(glm::vec3 p1, glm::vec3 p2, float radius, std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref) {
    glm::vec3 dir = p2 - p1;
    if (glm::length(dir) < 0.001f) return; 

    int base_index = vert_ref.size();
    int sides = 8; 

    glm::vec3 w = glm::normalize(dir);
    glm::vec3 u;
    if (std::abs(w.x) > 0.9f || std::abs(w.y) > 0.9f) {
        u = glm::normalize(glm::cross(w, {0.0f, 0.0f, 1.0f}));
    } else {
        u = glm::normalize(glm::cross(w, {1.0f, 0.0f, 0.0f}));
    }
    glm::vec3 v = glm::cross(w, u);

    for (int i = 0; i <= sides; ++i) {
        float angle = 2.0f * M_PI * i / sides;
        glm::vec3 offset = radius * (cosf(angle) * u + sinf(angle) * v);
        vert_ref.push_back(p1 + offset);
        vert_ref.push_back(p2 + offset);
    }

    for (int i = 0; i < sides; ++i) {
        int i0 = base_index + i * 2;
        int i1 = base_index + (i + 1) * 2;
        ind_ref.push_back(i0); ind_ref.push_back(i0 + 1); ind_ref.push_back(i1);
        ind_ref.push_back(i1); ind_ref.push_back(i0 + 1); ind_ref.push_back(i1 + 1);
    }
}

void generateHullAndDeck() {
    const int segments_per_curve = 30;
    const float bulwark_height = 0.8f;
    
    int num_points_per_ring = 2 * segments_per_curve - 1;

    for (size_t i = 0; i < blueprint_sections.size(); ++i) {
        for (int j = 0; j < segments_per_curve; ++j) {
            float t = (float)j / (segments_per_curve - 1);
            glm::vec2 p2d = blueprint_sections[i].getPoint(t);
            vertices.push_back(glm::vec3(p2d.x, p2d.y, blueprint_z_positions[i]));
        }
        for (int j = segments_per_curve - 2; j >= 0; --j) {
            float t = (float)j / (segments_per_curve - 1);
            glm::vec2 p2d = blueprint_sections[i].getPoint(t);
            vertices.push_back(glm::vec3(-p2d.x, p2d.y, blueprint_z_positions[i]));
        }
    }

    for (size_t i = 0; i < blueprint_sections.size() - 1; ++i) {
        int start1 = i * num_points_per_ring;
        int start2 = (i + 1) * num_points_per_ring;
        for (int j = 0; j < num_points_per_ring - 1; ++j) {
            indices.push_back(start1 + j); indices.push_back(start2 + j + 1); indices.push_back(start2 + j);
            indices.push_back(start1 + j); indices.push_back(start1 + j + 1); indices.push_back(start2 + j + 1);
        }
    }

    std::vector<int> top_edge_indices_right;
    std::vector<int> top_edge_indices_left;

    for (size_t i = 0; i < blueprint_sections.size(); ++i) {
        int ring_start_index = i * num_points_per_ring;
        top_edge_indices_right.push_back(ring_start_index + segments_per_curve - 1);
        top_edge_indices_left.push_back(ring_start_index + segments_per_curve);
    }

    for (size_t i = 0; i < blueprint_sections.size() -1; ++i) {
        indices.push_back(top_edge_indices_left[i]);
        indices.push_back(top_edge_indices_right[i+1]);
        indices.push_back(top_edge_indices_right[i]);
        indices.push_back(top_edge_indices_left[i]);
        indices.push_back(top_edge_indices_left[i+1]);
        indices.push_back(top_edge_indices_right[i+1]);
    }
    
    int bulwark_base_index_right = vertices.size();
    for(int idx : top_edge_indices_right) vertices.push_back(vertices[idx] + glm::vec3(0, bulwark_height, 0));
    int bulwark_base_index_left = vertices.size();
    for(int idx : top_edge_indices_left) vertices.push_back(vertices[idx] + glm::vec3(0, bulwark_height, 0));

    for (size_t i = 0; i < top_edge_indices_right.size() - 1; ++i) {
        int p1 = top_edge_indices_right[i]; int p2 = top_edge_indices_right[i+1];
        int p3 = bulwark_base_index_right + i; int p4 = bulwark_base_index_right + i + 1;
        indices.push_back(p1); indices.push_back(p2); indices.push_back(p4);
        indices.push_back(p1); indices.push_back(p4); indices.push_back(p3);
    }
     for (size_t i = 0; i < top_edge_indices_left.size() - 1; ++i) {
        int p1 = top_edge_indices_left[i]; int p2 = top_edge_indices_left[i+1];
        int p3 = bulwark_base_index_left + i; int p4 = bulwark_base_index_left + i + 1;
        indices.push_back(p1); indices.push_back(p4); indices.push_back(p2);
        indices.push_back(p1); indices.push_back(p3); indices.push_back(p4);
    }

    int stern_start_index = 0;
    for (int j = 0; j < num_points_per_ring - 1; ++j) {
        indices.push_back(stern_start_index + j);
        indices.push_back(stern_start_index);
        indices.push_back(stern_start_index + j + 1);
    }
    
    int bow_start_index = (blueprint_sections.size() - 1) * num_points_per_ring;
    glm::vec3 bow_tip = {0.0f, 5.0f, 30.0f};
    int bow_tip_index = vertices.size();
    vertices.push_back(bow_tip);
    for (int j = 0; j < num_points_per_ring - 1; ++j) {
        indices.push_back(bow_start_index + j);
        indices.push_back(bow_start_index + j + 1);
        indices.push_back(bow_tip_index);
    }
}

glm::vec3 generateBowspritPlatform() {
     float rope_radius = 0.05f;
    float post_radius = 0.06f;

    glm::vec3 base_pos = {0.0f, getDeckYAtZ(25.0f) + 0.2f, 25.0f}; 
    glm::vec3 tip_pos = {0.0f, base_pos.y + 1.0f, 45.0f};       

    generateRailSegment(base_pos, tip_pos, 0.3f, vertices, indices);

    float plank_width = 0.2f;
    float plank_gap = 0.05f;
    float platform_half_width = 1.0f;
    float plank_height = 0.08f;
    
    for (float z_offset = 1.0f; z_offset < 15.0f; z_offset += (plank_width + plank_gap)) {
        float t = z_offset / 15.0f;
        glm::vec3 plank_center = glm::mix(base_pos, tip_pos, t);
        float current_width = glm::mix(platform_half_width, 0.4f, t);

        generateBox({-current_width, plank_center.y - plank_height/2, plank_center.z}, {current_width - 0.15f, plank_height, plank_width}, vertices, indices);
        generateBox({0.15f, plank_center.y - plank_height/2, plank_center.z}, {current_width - 0.15f, plank_height, plank_width}, vertices, indices);
    }

    float rail_height = 1.0f;
    std::vector<glm::vec3> stanchion_tops;
    
    std::vector<glm::vec3> stanchion_bases;
    for (int i = 0; i < 5; ++i) {
        float t = i / 4.0f;
        float current_width = glm::mix(platform_half_width, 0.4f, t);
        glm::vec3 base_center = glm::mix(base_pos, tip_pos, t * 0.8f);
        stanchion_bases.push_back({-current_width, base_center.y, base_center.z});
    }
     for (int i = 4; i >= 0; --i) {
        float t = i / 4.0f;
        float current_width = glm::mix(platform_half_width, 0.4f, t);
        glm::vec3 base_center = glm::mix(base_pos, tip_pos, t * 0.8f);
        stanchion_bases.push_back({current_width, base_center.y, base_center.z});
    }

    for(const auto& base : stanchion_bases) {
        generateCylinder(base, rail_height, post_radius, 6, vertices, indices);
        stanchion_tops.push_back(base + glm::vec3(0, rail_height, 0));
    }
    
    for(size_t i = 0; i < stanchion_tops.size() - 1; ++i) {
        glm::vec3 current_top = stanchion_tops[i];
        glm::vec3 next_top = stanchion_tops[i+1];
        
        generateRailSegment(current_top, next_top, rope_radius, vertices, indices); // Pagar atas
        generateRailSegment(current_top - glm::vec3(0, rail_height / 2, 0), next_top - glm::vec3(0, rail_height / 2, 0), rope_radius, vertices, indices); // Pagar tengah
    }

    return tip_pos; 
}

//GENERATE SPARS UNTUK RIG PHINISI (DINAMIS)
MastGeometry generateSingleMast(glm::vec3 base_pos, float height, bool is_foremast) {
    MastGeometry geo;
    geo.base = base_pos;
    geo.top = base_pos + glm::vec3(0, height, 0);

    float boom_radius = 0.25f;
    float gaff_radius = 0.22f;

    if (is_foremast) { 
        geo.boom_base = geo.base + glm::vec3(0, 4.5f, 0);
        geo.boom_tip  = {0.0f, geo.boom_base.y + 1.0f, -8.0f};
        geo.gaff_base = geo.base + glm::vec3(0, height * 0.88f, 0);
        geo.gaff_tip  = {0.0f, geo.gaff_base.y + 6.0f, -4.0f};
        geo.jib_boom_base = geo.base + glm::vec3(0, 2.0f, 0);
        geo.jib_boom_tip  = {0.0f, geo.jib_boom_base.y + 0.2f, 26.0f};
    } else {
        geo.boom_base = geo.base + glm::vec3(0, 4.0f, 0);
        geo.boom_tip  = {0.0f, geo.boom_base.y + 1.0f, -28.0f};
        geo.gaff_base = geo.base + glm::vec3(0, height * 0.85f, 0);
        geo.gaff_tip  = {0.0f, geo.gaff_base.y + 5.0f, -25.0f};
    }
    
    generateCylinder(geo.base, height, 0.45f, 16, vertices, indices);
    generateRailSegment(geo.boom_base, geo.boom_tip, boom_radius, vertices, indices);
    generateRailSegment(geo.gaff_base, geo.gaff_tip, gaff_radius, vertices, indices);
    generateRailSegment(geo.jib_boom_base, geo.jib_boom_tip, boom_radius * 0.8f, vertices, indices);
    
    return geo;
}



void generateSails() {
    glm::vec3 foremast_pos = {0.0f, getDeckYAtZ(12.0f), 12.0f};
    float foremast_height = 45.0f;
    glm::vec3 foremast_top = foremast_pos + glm::vec3(0, foremast_height, 0);

    glm::vec3 mainmast_pos = {0.0f, getDeckYAtZ(-12.0f), -12.0f};
    float mainmast_height = 40.0f;
    glm::vec3 mainmast_top = mainmast_pos + glm::vec3(0, mainmast_height, 0);

    glm::vec3 bowsprit_tip = {0.0f, getDeckYAtZ(25.0f) + 2.8f, 45.0f}; 

    glm::vec3 mainmast_boom_base = mainmast_pos + glm::vec3(0, 4.0f, 0);
    glm::vec3 mainmast_boom_tip  = {0.0f, mainmast_boom_base.y + 1.0f, -24.0f};
    glm::vec3 mainmast_gaff_base = mainmast_pos + glm::vec3(0, mainmast_height * 0.85f, 0);
    glm::vec3 mainmast_gaff_tip  = {0.0f, mainmast_gaff_base.y + 5.0f, -20.0f};

    glm::vec3 foremast_boom_base = foremast_pos + glm::vec3(0, 4.5f, 0);
    glm::vec3 foremast_boom_tip  = {0.0f, foremast_boom_base.y + 1.0f, -3.0f};
    glm::vec3 foremast_gaff_base = foremast_pos + glm::vec3(0, foremast_height * 0.88f, 0);
    glm::vec3 foremast_gaff_tip  = {0.0f, foremast_gaff_base.y + 6.0f, 2.0f};
    
    int divisions = 20; 

    generateSailMesh(mainmast_gaff_tip, mainmast_boom_tip, mainmast_boom_base, mainmast_gaff_base, divisions, sail_vertices, sail_indices);

    generateSailMesh(foremast_gaff_tip, foremast_boom_tip, foremast_boom_base, foremast_gaff_base, divisions, sail_vertices, sail_indices);

    generateSailMesh(mainmast_top, mainmast_gaff_tip, mainmast_gaff_base, mainmast_top, divisions, sail_vertices, sail_indices);

    generateSailMesh(foremast_top, foremast_gaff_tip, foremast_gaff_base, foremast_top, divisions, sail_vertices, sail_indices);

    glm::vec3 jib_clew = {0.0f, getDeckYAtZ(20.0f) + 1.0f, 20.0f};
    
    generateSailMesh(foremast_top, jib_clew, bowsprit_tip, bowsprit_tip, 12, sail_vertices, sail_indices);
    generateSailMesh({foremast_top.x, foremast_top.y - 15.0f, foremast_pos.z}, jib_clew, {bowsprit_tip.x, bowsprit_tip.y-1.0f, bowsprit_tip.z - 12.0f}, {bowsprit_tip.x, bowsprit_tip.y-1.0f, bowsprit_tip.z - 12.0f}, 12, sail_vertices, sail_indices);
}

// GENERATE TANGGA TIANG SECARA DINAMIS
void generateMastLadder(glm::vec3 mast_base, float mast_height, float mast_radius, 
                      const glm::vec3& orientation, 
                      std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref) {
    std::cout << "    -> Membangun tangga kayu untuk tiang dengan orientasi kustom..." << std::endl;

    const float LADDER_WIDTH = 0.5f;
    const float RAIL_RADIUS = 0.04f;
    const float RUNG_RADIUS = 0.03f;
    const float RUNG_SPACING = 0.4f;
    const float STANDOFF_DISTANCE = 0.1f;

    const glm::vec3 standoff_vec = glm::normalize(orientation); 
    const glm::vec3 side_vec = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), standoff_vec));

    glm::vec3 ladder_center_base = mast_base + standoff_vec * (mast_radius + STANDOFF_DISTANCE);
    glm::vec3 left_rail_base = ladder_center_base - side_vec * (LADDER_WIDTH / 2.0f);
    glm::vec3 right_rail_base = ladder_center_base + side_vec * (LADDER_WIDTH / 2.0f);

    glm::vec3 left_rail_top = left_rail_base + glm::vec3(0, mast_height, 0);
    glm::vec3 right_rail_top = right_rail_base + glm::vec3(0, mast_height, 0);

    generateRailSegment(left_rail_base, left_rail_top, RAIL_RADIUS, vert_ref, ind_ref);
    generateRailSegment(right_rail_base, right_rail_top, RAIL_RADIUS, vert_ref, ind_ref);

    int num_rungs = static_cast<int>(mast_height / RUNG_SPACING);
    for (int i = 1; i < num_rungs; ++i) {
        float current_height = i * RUNG_SPACING;
        glm::vec3 rung_left_pos = left_rail_base + glm::vec3(0, current_height, 0);
        glm::vec3 rung_right_pos = right_rail_base + glm::vec3(0, current_height, 0);
        
        generateRailSegment(rung_left_pos, rung_right_pos, RUNG_RADIUS, vert_ref, ind_ref);
    }
}

// GENERATE JARING BOWSPRIT DENGAN EFEK LENDUTAN
void generateSideNetting(bool is_right_side, std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref) {
    std::cout << "    -> Merajut jaring samping yang realistis..." << std::endl;

    const float NET_LENGTH = 15.0f;     
    const float NET_HEIGHT = 4.0f;     
    const int NUM_ROPES_Z = 20;         
    const int NUM_ROPES_Y = 8;          
    const float ROPE_RADIUS = 0.03f;
    const float SAG_Y = 0.5f;           
    const float SAG_X = 0.4f;           
    const float RANDOMNESS = 0.05f;     
    
    const float start_z = -5.0f;
    float side_multiplier = is_right_side ? 1.0f : -1.0f;

    auto getNetNode = [&](float t_z, float t_y) {

        float current_z = start_z + t_z * NET_LENGTH;
        
        glm::vec3 anchor_point;
        anchor_point.z = current_z;
        anchor_point.y = getDeckYAtZ(current_z) + 0.8f; 
        anchor_point.x = getDeckWidthAtZ(current_z) * side_multiplier;

        glm::vec3 point = anchor_point - glm::vec3(0.0f, t_y * NET_HEIGHT, 0.0f);

        point.y -= SAG_Y * sinf(t_z * M_PI);

        float outward_sag = SAG_X * sinf(t_y * M_PI / 2.0f);
        point.x += outward_sag * side_multiplier;

        if (RANDOMNESS > 0.0f) {
            auto random_offset = [&]() { return (static_cast<float>(rand()) / RAND_MAX - 0.5f) * 2.0f * RANDOMNESS; };
            point += glm::vec3(random_offset(), random_offset(), random_offset());
        }

        return point;
    };

    for (int i = 0; i <= NUM_ROPES_Y; ++i) {
        float t_y = static_cast<float>(i) / NUM_ROPES_Y;
        for (int j = 0; j < NUM_ROPES_Z; ++j) {
            float t_z1 = static_cast<float>(j) / NUM_ROPES_Z;
            float t_z2 = static_cast<float>(j + 1) / NUM_ROPES_Z;
            generateRailSegment(getNetNode(t_z1, t_y), getNetNode(t_z2, t_y), ROPE_RADIUS, vert_ref, ind_ref);
        }
    }
    for (int i = 0; i <= NUM_ROPES_Z; ++i) {
        float t_z = static_cast<float>(i) / NUM_ROPES_Z;
        for (int j = 0; j < NUM_ROPES_Y; ++j) {
            float t_y1 = static_cast<float>(j) / NUM_ROPES_Y;
            float t_y2 = static_cast<float>(j + 1) / NUM_ROPES_Y;
            generateRailSegment(getNetNode(t_z, t_y1), getNetNode(t_z, t_y2), ROPE_RADIUS, vert_ref, ind_ref);
        }
    }
}

// GENERATE SATU BUAH TONG KAYU DETAIL
void generateBarrel(glm::vec3 base_pos, float height, float mid_radius, float end_radius, int segments, int slices, std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref)
{
    glm::vec2 P0(end_radius, 0.0f);
    glm::vec2 P1(mid_radius, height * 0.25f);
    glm::vec2 P2(mid_radius, height * 0.75f);
    glm::vec2 P3(end_radius, height);

    auto getProfilePoint = [&](float t) {
        float u = 1.0f - t;
        float tt = t * t;
        float uu = u * u;
        float uuu = uu * u;
        float ttt = tt * t;
        return uuu * P0 + 3.0f * uu * t * P1 + 3.0f * u * tt * P2 + ttt * P3;
    };

    int base_index = vert_ref.size();
    for (int i = 0; i <= slices; ++i) {
        float t = static_cast<float>(i) / slices;
        glm::vec2 profile_pt = getProfilePoint(t);
        for (int j = 0; j <= segments; ++j) {
            float angle = 2.0f * M_PI * j / segments;
            glm::vec3 point(profile_pt.x * cosf(angle), profile_pt.y, profile_pt.x * sinf(angle));
            vert_ref.push_back(base_pos + point);
        }
    }

    for (int i = 0; i < slices; ++i) {
        for (int j = 0; j < segments; ++j) {
            int row1 = i * (segments + 1);
            int row2 = (i + 1) * (segments + 1);
            ind_ref.push_back(base_index + row1 + j);
            ind_ref.push_back(base_index + row2 + j + 1);
            ind_ref.push_back(base_index + row1 + j + 1);
            
            ind_ref.push_back(base_index + row1 + j);
            ind_ref.push_back(base_index + row2 + j);
            ind_ref.push_back(base_index + row2 + j + 1);
        }
    }

    int top_center_idx = vert_ref.size();
    vert_ref.push_back(base_pos + glm::vec3(0, height, 0));
    int top_ring_start = base_index + slices * (segments + 1);
    for (int i = 0; i < segments; i++) {
        ind_ref.push_back(top_center_idx);
        ind_ref.push_back(top_ring_start + i + 1);
        ind_ref.push_back(top_ring_start + i);
    }
    
    float hoop_heights[] = {height * 0.05f, height * 0.15f, height * 0.85f, height * 0.95f};
    float hoop_thickness = 0.05f;
    for(float h : hoop_heights) {
        float t = h / height;
        float r = getProfilePoint(t).x + 0.01f;
        generateCylinder(base_pos + glm::vec3(0, h - (hoop_thickness/2), 0), hoop_thickness, r, segments, vertices, indices);
    }
}

// GENERATE AREA MAKAN DI ATAP DECKHOUSE
void generateRooftopDining() {
    std::cout << "    -> Membangun area makan mewah di atap..." << std::endl;

    glm::vec3 deckhouse_base_pos = {-3.0f, getDeckYAtZ(2.5f), -4.0f};
    glm::vec3 deckhouse_body_size = {6.0f, 2.5f, 9.0f};
    glm::vec3 roof_pos = deckhouse_base_pos + glm::vec3(-0.2f, deckhouse_body_size.y, -0.2f);
    glm::vec3 roof_size = {deckhouse_body_size.x + 0.4f, 0.2f, deckhouse_body_size.z + 0.4f};
    float rooftop_y = roof_pos.y + roof_size.y;

    float table_height = 0.8f;
    float table_length = 4.0f;
    float table_width = 1.5f;
    glm::vec3 table_center_pos = {roof_pos.x + roof_size.x / 2.0f, rooftop_y, roof_pos.z + roof_size.z / 2.0f};
    glm::vec3 tabletop_pos = table_center_pos + glm::vec3(-table_length/2, table_height, -table_width/2);
    generateBox(tabletop_pos, {table_length, 0.1f, table_width}, vertices, indices);
    generateBox(tabletop_pos + glm::vec3(0.1f, -table_height, 0.1f), {0.2f, table_height, 0.2f}, vertices, indices);
    generateBox(tabletop_pos + glm::vec3(table_length - 0.3f, -table_height, 0.1f), {0.2f, table_height, 0.2f}, vertices, indices);
    generateBox(tabletop_pos + glm::vec3(0.1f, -table_height, table_width - 0.3f), {0.2f, table_height, 0.2f}, vertices, indices);
    generateBox(tabletop_pos + glm::vec3(table_length - 0.3f, -table_height, table_width - 0.3f), {0.2f, table_height, 0.2f}, vertices, indices);

    float bench_height = 0.5f;
    glm::vec3 bench1_pos = tabletop_pos + glm::vec3(0, -table_height, -table_width/2 - 0.5f);
    generateBox(bench1_pos, {table_length, 0.1f, 0.4f}, vertices, indices);
    generateBox(bench1_pos + glm::vec3(0,0,0), {table_length, -bench_height, 0.1f}, vertices, indices);
    generateBox(bench1_pos + glm::vec3(0,0,0.3f), {table_length, -bench_height, 0.1f}, vertices, indices);
    glm::vec3 bench2_pos = tabletop_pos + glm::vec3(0, -table_height, table_width + 0.1f);
    generateBox(bench2_pos, {table_length, 0.1f, 0.4f}, vertices, indices);
    generateBox(bench2_pos + glm::vec3(0,0,0), {table_length, -bench_height, 0.1f}, vertices, indices);
    generateBox(bench2_pos + glm::vec3(0,0,0.3f), {table_length, -bench_height, 0.1f}, vertices, indices);

    float pergola_height = 1.5f; 
    
    glm::vec3 p1 = roof_pos + glm::vec3(0.5f, roof_size.y, 0.5f);
    glm::vec3 p2 = roof_pos + glm::vec3(roof_size.x - 0.5f, roof_size.y, 0.5f);
    glm::vec3 p3 = roof_pos + glm::vec3(roof_size.x - 0.5f, roof_size.y, roof_size.z - 0.5f);
    glm::vec3 p4 = roof_pos + glm::vec3(0.5f, roof_size.y, roof_size.z - 0.5f);

    float cover_height = pergola_height + 0.1f; 
    generateBox(roof_pos + glm::vec3(0.3f, roof_size.y + pergola_height - 0.1f, 0.3f),
               {roof_size.x - 0.6f, 0.1f, roof_size.z - 0.6f}, 
               vertices, indices);

    float cross_height = pergola_height - 0.2f;
    for(float z = 0.5f; z < roof_size.z - 0.5f; z += 1.0f) {
        generateBox(roof_pos + glm::vec3(0.3f, roof_size.y + cross_height, z),
                   {roof_size.x - 0.6f, 0.05f, 0.05f},
                   vertices, indices);
    }
    
    generateCylinder(p1, pergola_height, 0.08f, 8, vertices, indices);
    generateCylinder(p2, pergola_height, 0.08f, 8, vertices, indices);
    generateCylinder(p3, pergola_height, 0.08f, 8, vertices, indices);
    generateCylinder(p4, pergola_height, 0.08f, 8, vertices, indices);

    glm::vec3 top1 = p1 + glm::vec3(0,pergola_height,0);
    glm::vec3 top2 = p2 + glm::vec3(0,pergola_height,0);
    glm::vec3 top3 = p3 + glm::vec3(0,pergola_height,0);
    glm::vec3 top4 = p4 + glm::vec3(0,pergola_height,0);
    generateBox(top1 + glm::vec3(0,0,-0.1f), {glm::distance(top1,top2), 0.1f, 0.2f}, vertices, indices);
    generateBox(top4 + glm::vec3(0,0,-0.1f), {glm::distance(top4,top3), 0.1f, 0.2f}, vertices, indices);
    generateBox(top1 + glm::vec3(-0.1f,0,0), {0.2f, 0.1f, glm::distance(top1,top4)}, vertices, indices);
    generateBox(top2 + glm::vec3(-0.1f,0,0), {0.2f, 0.1f, glm::distance(top2,top3)}, vertices, indices);
}

// MENGATUR POSISI Tumpukan TONG
void generateBarrelCluster(glm::vec3 cluster_base_pos) { 
    std::cout << "    -> Menambahkan tumpukan tong di dek..." << std::endl;
    
    float h = 1.2f, mid_r = 0.7f, end_r = 0.6f;
    int seg = 20, sli = 10;

    generateBarrel(cluster_base_pos, h, mid_r, end_r, seg, sli, vertices, indices);
    
    generateBarrel(cluster_base_pos + glm::vec3(1.0f, 0, 0.5f), h, mid_r, end_r, seg, sli, vertices, indices);

    int pre_tilt_v_count = vertices.size();
    glm::vec3 tilted_barrel_pos = cluster_base_pos + glm::vec3(-1.0f, 0, 0.5f); 
    generateBarrel(tilted_barrel_pos, h, mid_r, end_r, seg, sli, vertices, indices);
    int post_tilt_v_count = vertices.size();

    glm::mat4 rotation_matrix = glm::rotate(glm::mat4(1.0f), glm::radians(15.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    for (int i = pre_tilt_v_count; i < post_tilt_v_count; ++i) {
        vertices[i] -= tilted_barrel_pos;
        vertices[i] = glm::vec3(rotation_matrix * glm::vec4(vertices[i], 1.0f));
        vertices[i] += tilted_barrel_pos; 
    }
}

// GENERATE GULUNGAN TALI SPIRAL
void generateCoiledRope(glm::vec3 center_pos, float max_radius, int coils, float rope_radius, std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref) {
    std::cout << "    -> Menggulung tali di dek..." << std::endl;
    
    int segments_per_coil = 30;
    int total_segments = coils * segments_per_coil;
    float angle_step = 2.0f * M_PI / segments_per_coil;
    
    glm::vec3 p1 = center_pos;

    for (int i = 1; i <= total_segments; ++i) {
        float current_angle = i * angle_step;
        float current_radius = (static_cast<float>(i) / total_segments) * max_radius;
        
        glm::vec3 p2 = center_pos + glm::vec3(current_radius * cosf(current_angle), i * 0.001f, current_radius * sinf(current_angle));
        
        glm::vec3 dir = p2 - p1;
        if (glm::length(dir) > 0.001f) {
            generateCylinder(p1, glm::distance(p1, p2), rope_radius, 4, vert_ref, ind_ref, glm::normalize(dir));
        }
        
        p1 = p2;
    }
}

// MENEMPATKAN BEBERAPA GULUNGAN TALI
void placeCoiledRopes() {
    glm::vec3 mainmast_base = {0.0f, getDeckYAtZ(-12.0f), -12.0f};
    generateCoiledRope(mainmast_base + glm::vec3(1.5f, 0.0f, 0.0f), 0.8f, 5, 0.04f, vertices, indices);

    glm::vec3 foremast_base = {0.0f, getDeckYAtZ(12.0f), 12.0f};
    generateCoiledRope(foremast_base + glm::vec3(-1.5f, 0.0f, 0.0f), 0.8f, 5, 0.04f, vertices, indices);
}

void generateTorus(glm::vec3 center, float major_radius, float minor_radius, int major_segments, int minor_segments, std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref) {
    std::cout << "    -> Membuat pelampung penyelamat..." << std::endl;
    int base_index = vert_ref.size();

    // Generate vertices
    for (int i = 0; i <= major_segments; i++) {
        float theta = 2.0f * M_PI * static_cast<float>(i) / major_segments;
        glm::vec3 major_pos(major_radius * cosf(theta), 0, major_radius * sinf(theta));

        for (int j = 0; j <= minor_segments; j++) {
            float phi = 2.0f * M_PI * static_cast<float>(j) / minor_segments;
            glm::vec3 normal = glm::normalize(glm::vec3(cosf(theta) * cosf(phi), sinf(phi), sinf(theta) * cosf(phi)));
            vert_ref.push_back(center + major_pos + normal * minor_radius);
        }
    }

    for (int i = 0; i < major_segments; i++) {
        for (int j = 0; j < minor_segments; j++) {
            int p1 = i * (minor_segments + 1) + j;
            int p2 = i * (minor_segments + 1) + (j + 1);
            int p3 = (i + 1) * (minor_segments + 1) + j;
            int p4 = (i + 1) * (minor_segments + 1) + (j + 1);
            ind_ref.push_back(base_index + p1); ind_ref.push_back(base_index + p3); ind_ref.push_back(base_index + p2);
            ind_ref.push_back(base_index + p2); ind_ref.push_back(base_index + p3); ind_ref.push_back(base_index + p4);
        }
    }
}

// MENEMPATKAN PELAMPUNG PENYELAMAT
void placeLifebuoys() {
    float deck_y = getDeckYAtZ(-15.0f) + 0.8f;
    float deck_width = getDeckWidthAtZ(-15.0f);

    glm::vec3 pos1 = {-deck_width * 0.9f - 0.1f, deck_y + 1.5f, -17.0f};
    
    int pre_rot_v_count = vertices.size();
    generateTorus(pos1, 0.4f, 0.08f, 30, 15, vertices, indices);
    int post_rot_v_count = vertices.size();

    glm::mat4 rotation_matrix = glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    for (int i = pre_rot_v_count; i < post_rot_v_count; ++i) {
        vertices[i] -= pos1;
        vertices[i] = glm::vec3(rotation_matrix * glm::vec4(vertices[i], 1.0f));
        vertices[i] += pos1;
    }
}
// GENERATE SILINDER DENGAN ORIENTASI APAPUN
void generateOrientedCylinder(glm::vec3 start_pos, glm::vec3 end_pos, float radius, int sides, std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref) {
    glm::vec3 dir = end_pos - start_pos;
    if (glm::length(dir) < 0.001f) return;

    int base_index = vert_ref.size();
    
    glm::vec3 w = glm::normalize(dir);
    glm::vec3 u = (std::abs(w.y) > 0.99f) ? glm::normalize(glm::cross(w, {0.0f, 0.0f, 1.0f})) : glm::normalize(glm::cross(w, {0.0f, 1.0f, 0.0f}));
    glm::vec3 v = glm::cross(w, u);

    for (int i = 0; i <= sides; ++i) {
        float angle = 2.0f * M_PI * i / sides;
        glm::vec3 offset = radius * (cosf(angle) * u + sinf(angle) * v);
        vert_ref.push_back(start_pos + offset);
    }
    for (int i = 0; i <= sides; ++i) {
        float angle = 2.0f * M_PI * i / sides;
        glm::vec3 offset = radius * (cosf(angle) * u + sinf(angle) * v);
        vert_ref.push_back(end_pos + offset);
    }

    int ring1_start = base_index;
    int ring2_start = base_index + (sides + 1);
    for (int i = 0; i < sides; ++i) {
        ind_ref.push_back(ring1_start + i); ind_ref.push_back(ring2_start + i); ind_ref.push_back(ring1_start + i + 1);
        ind_ref.push_back(ring2_start + i); ind_ref.push_back(ring2_start + i + 1); ind_ref.push_back(ring1_start + i + 1);
    }

    int start_center_idx = vert_ref.size(); vert_ref.push_back(start_pos);
    int end_center_idx = vert_ref.size(); vert_ref.push_back(end_pos);
    for (int i = 0; i < sides; i++) {
        ind_ref.push_back(start_center_idx); ind_ref.push_back(ring1_start + i); ind_ref.push_back(ring1_start + i + 1);
        ind_ref.push_back(end_center_idx); ind_ref.push_back(ring2_start + i + 1); ind_ref.push_back(ring2_start + i);
    }
}
//  GENERATE DECKHOUSE MEWAH DAN DETAIL
void generateDeckhouse() {
    std::cout << "    -> Finalisasi Deckhouse: Memperbaiki jendela dan menambah tangga..." << std::endl;

    glm::vec3 base_pos = {-3.0f, getDeckYAtZ(2.5f), -4.0f};
    glm::vec3 body_size = {6.0f, 2.5f, 9.0f};
    float wall_thickness = 0.15f;

    generateBox(base_pos, {wall_thickness, body_size.y, body_size.z}, vertices, indices);
    generateBox(base_pos + glm::vec3(wall_thickness, 0, 0), {body_size.x - wall_thickness*2, body_size.y, wall_thickness}, vertices, indices);
    generateBox(base_pos + glm::vec3(wall_thickness, 0, body_size.z - wall_thickness), {body_size.x - wall_thickness*2, body_size.y, wall_thickness}, vertices, indices);
    
    float door_width = 1.2f;
    float door_height = 2.0f;
    float door_pos_z = 3.5f;
    glm::vec3 right_wall_pos = base_pos + glm::vec3(body_size.x - wall_thickness, 0, 0);
    generateBox(right_wall_pos, {wall_thickness, body_size.y, door_pos_z}, vertices, indices);
    generateBox(right_wall_pos + glm::vec3(0,0,door_pos_z + door_width), {wall_thickness, body_size.y, body_size.z - (door_pos_z + door_width)}, vertices, indices);
    generateBox(right_wall_pos + glm::vec3(0, door_height, door_pos_z), {wall_thickness, body_size.y - door_height, door_width}, vertices, indices);

    glm::vec3 door_pos = right_wall_pos + glm::vec3(0, 0, door_pos_z);
    generateBox(door_pos, {wall_thickness, door_height, door_width}, vertices, indices);
    glm::vec3 handle_pos = door_pos + glm::vec3(wall_thickness/2, door_height/2, door_width - 0.2f);
    generateRailSegment(handle_pos, handle_pos + glm::vec3(0.1f, 0, 0), 0.03f, vertices, indices);

    float porthole_radius = 0.4f;
    for(int i=0; i<2; ++i) {
        glm::vec3 p_center = base_pos + glm::vec3(wall_thickness/2, 1.5f, 2.5f + i * 4.0f);
        generateOrientedCylinder(p_center - glm::vec3(wall_thickness,0,0), p_center + glm::vec3(wall_thickness,0,0), porthole_radius, 16, vertices, indices);
    }
    glm::vec3 p_center_right = right_wall_pos + glm::vec3(-wall_thickness/2, 1.5f, body_size.z - 2.5f);
    generateOrientedCylinder(p_center_right - glm::vec3(wall_thickness,0,0), p_center_right + glm::vec3(wall_thickness,0,0), porthole_radius, 16, vertices, indices);


    glm::vec3 roof_pos = base_pos + glm::vec3(-0.2f, body_size.y, -0.2f);
    glm::vec3 roof_size = {body_size.x + 0.4f, 0.2f, body_size.z + 0.4f};
    generateBox(roof_pos, roof_size, vertices, indices); 
    float rail_height = 0.5f; float post_radius = 0.04f;
    std::vector<glm::vec3> rail_posts;
    rail_posts.push_back(roof_pos + glm::vec3(post_radius, roof_size.y, post_radius));
    rail_posts.push_back(roof_pos + glm::vec3(roof_size.x - post_radius, roof_size.y, post_radius));
    rail_posts.push_back(roof_pos + glm::vec3(roof_size.x - post_radius, roof_size.y, roof_size.z - post_radius));
    rail_posts.push_back(roof_pos + glm::vec3(post_radius, roof_size.y, roof_size.z - post_radius));
    for(const auto& post_base : rail_posts) { generateCylinder(post_base, rail_height, post_radius, 8, vertices, indices); }
    for(size_t i=0; i<rail_posts.size(); ++i) {
        glm::vec3 p1 = rail_posts[i] + glm::vec3(0,rail_height,0);
        glm::vec3 p2 = rail_posts[(i+1) % rail_posts.size()] + glm::vec3(0,rail_height,0);
        generateRailSegment(p1, p2, post_radius * 0.8f, vertices, indices);
    }
    
    float ladder_width = 0.6f;
    float rung_radius = 0.03f;
    glm::vec3 ladder_back_base = base_pos + glm::vec3(body_size.x/2, 0, -0.1f);
    glm::vec3 ladder_back_top = ladder_back_base + glm::vec3(0, body_size.y, 0);
    generateRailSegment(ladder_back_base - glm::vec3(ladder_width/2,0,0), ladder_back_top - glm::vec3(ladder_width/2,0,0), rung_radius, vertices, indices);
    generateRailSegment(ladder_back_base + glm::vec3(ladder_width/2,0,0), ladder_back_top + glm::vec3(ladder_width/2,0,0), rung_radius, vertices, indices);
    for(float y=0.3f; y < body_size.y; y += 0.4f) {
        generateRailSegment(ladder_back_base - glm::vec3(ladder_width/2,-y,0), ladder_back_base + glm::vec3(ladder_width/2,y,0), rung_radius, vertices, indices);
    }
    glm::vec3 ladder_front_base = base_pos + glm::vec3(body_size.x/2, 0, body_size.z + 0.1f);
    glm::vec3 ladder_front_top = ladder_front_base + glm::vec3(0, body_size.y, 0);
    generateRailSegment(ladder_front_base - glm::vec3(ladder_width/2,0,0), ladder_front_top - glm::vec3(ladder_width/2,0,0), rung_radius, vertices, indices);
    generateRailSegment(ladder_front_base + glm::vec3(ladder_width/2,0,0), ladder_front_top + glm::vec3(ladder_width/2,0,0), rung_radius, vertices, indices);
    for(float y=0.3f; y < body_size.y; y += 0.4f) {
        generateRailSegment(ladder_front_base - glm::vec3(ladder_width/2,-y,0), ladder_front_base + glm::vec3(ladder_width/2,y,0), rung_radius, vertices, indices);
    }
}

// GENERATE RIGGING PHINISI (PRESISI & DINAMIS)
void generateRigging_Phinisi(const MastGeometry& foremast, const MastGeometry& mainmast, const glm::vec3& bowsprit_tip) {
    float rope_radius = 0.05f;
    float thin_rope_radius = 0.03f;
    
    std::cout << "  - Memasang Standing Rigging..." << std::endl;
    generateRailSegment(foremast.top, bowsprit_tip, rope_radius, vertices, indices);
    generateRailSegment(mainmast.top, foremast.top, rope_radius, vertices, indices);
    generateRailSegment(mainmast.top, {0.0f, getDeckYAtZ(-25.0f) + 2.0f, -25.0f}, rope_radius, vertices, indices);

    auto generateShrouds = [&](const MastGeometry& mast, float deck_width){
        float hull_attachment_width = deck_width + 0.3f;
        float hull_attachment_y = mast.base.y - 1.0f;

        glm::vec3 deck_left = {-hull_attachment_width, hull_attachment_y, mast.base.z};
        glm::vec3 deck_right = {hull_attachment_width, hull_attachment_y, mast.base.z};
        
        generateBox({deck_left.x - 0.1f, deck_left.y - 0.2f, deck_left.z - 0.5f}, {0.2f, 0.4f, 1.0f}, vertices, indices);
        generateBox({deck_right.x - 0.1f, deck_right.y - 0.2f, deck_right.z - 0.5f}, {0.2f, 0.4f, 1.0f}, vertices, indices);
        
        generateRailSegment(mast.top, deck_left, thin_rope_radius, vertices, indices);
        generateRailSegment(mast.top, deck_right, thin_rope_radius, vertices, indices);
    };
    generateShrouds(foremast, getDeckWidthAtZ(foremast.base.z));
    generateShrouds(mainmast, getDeckWidthAtZ(mainmast.base.z));


    std::cout << "  - Memasang Running Rigging..." << std::endl;
    generateRailSegment(foremast.gaff_tip, foremast.top, thin_rope_radius, vertices, indices);
    generateRailSegment(mainmast.gaff_tip, mainmast.top, thin_rope_radius, vertices, indices);

    generateRailSegment(foremast.gaff_base, foremast.top, thin_rope_radius, vertices, indices);
    generateRailSegment(mainmast.gaff_base, mainmast.top, thin_rope_radius, vertices, indices);

    glm::vec3 sheet_point_deck_left = {-2.0f, getDeckYAtZ(-22.0f)+0.8f, -22.0f};
    glm::vec3 sheet_point_deck_right = {2.0f, getDeckYAtZ(-22.0f)+0.8f, -22.0f};
    generateRailSegment(mainmast.boom_tip, sheet_point_deck_left, thin_rope_radius, vertices, indices);
    generateRailSegment(mainmast.boom_tip, sheet_point_deck_right, thin_rope_radius, vertices, indices);
    
    generateRailSegment(foremast.boom_tip, foremast.top, thin_rope_radius, vertices, indices);
    generateRailSegment(mainmast.boom_tip, mainmast.top, thin_rope_radius, vertices, indices);
}

void generateRudder() {
    glm::vec3 rudder_top_pos = {0.0f, 7.0f, -22.5f};
    glm::vec3 rudder_size = {0.4f, 6.0f, 2.5f};
    
    int base_index = vertices.size();
    vertices.push_back(rudder_top_pos + glm::vec3(-rudder_size.x/2, 0, 0));
    vertices.push_back(rudder_top_pos + glm::vec3( rudder_size.x/2, 0, 0));
    vertices.push_back(rudder_top_pos + glm::vec3( rudder_size.x/2, -rudder_size.y, -rudder_size.z * 0.2f));
    vertices.push_back(rudder_top_pos + glm::vec3(-rudder_size.x/2, -rudder_size.y, -rudder_size.z * 0.2f));
    vertices.push_back(rudder_top_pos + glm::vec3(-rudder_size.x/2, 0, -rudder_size.z));
    vertices.push_back(rudder_top_pos + glm::vec3( rudder_size.x/2, 0, -rudder_size.z));
    vertices.push_back(rudder_top_pos + glm::vec3( rudder_size.x/2, -rudder_size.y * 0.8f, -rudder_size.z));
    vertices.push_back(rudder_top_pos + glm::vec3(-rudder_size.x/2, -rudder_size.y * 0.8f, -rudder_size.z));

    int indices_data[] = {0,1,2, 0,2,3, 4,5,6, 4,6,7, 0,4,7, 0,7,3, 1,5,6, 1,6,2, 0,1,5, 0,5,4, 3,2,6, 3,6,7};
    for(int i=0; i<36; ++i) indices.push_back(base_index + indices_data[i]);
}

void generateDeckRailings() {
    float rail_height = 0.8f;
    float post_radius = 0.05f;
    std::vector<glm::vec3> right_posts, left_posts;
    
    for(float z = -20.0f; z < 25.0f; z += 2.5f) {
        float deck_y = getDeckYAtZ(z) + 0.8f;
        float x_pos = getDeckWidthAtZ(z);

        glm::vec3 right_post_base = {x_pos, deck_y, z};
        glm::vec3 left_post_base = {-x_pos, deck_y, z};
        generateCylinder(right_post_base, rail_height, post_radius, 6, vertices, indices);
        generateCylinder(left_post_base, rail_height, post_radius, 6, vertices, indices);
        right_posts.push_back(right_post_base + glm::vec3(0, rail_height, 0));
        left_posts.push_back(left_post_base + glm::vec3(0, rail_height, 0));
    }
    for(size_t i = 0; i < right_posts.size() - 1; ++i) {
        generateRailSegment(right_posts[i], right_posts[i+1], post_radius, vertices, indices);
        generateRailSegment(left_posts[i], left_posts[i+1], post_radius, vertices, indices);
    }
    generateRailSegment(right_posts.back(), left_posts.back(), post_radius, vertices, indices);
}

// GENERATE PAPAN DEK DINAMIS SESUAI BENTUK LAMBUNG
void generateDeckPlanking() {
    std::cout << "    -> Memasang papan dek secara dinamis..." << std::endl;
    
    const float plank_width = 0.2f;
    const float plank_height = 0.05f;
    const float plank_segment_length = 1.0f;

    for (float z = -22.0f; z < 25.0f; z += plank_segment_length) {
        
        float current_deck_width = getDeckWidthAtZ(z);
        float current_deck_y = getDeckYAtZ(z) - plank_height;

        for (float x = -current_deck_width; x < current_deck_width; x += plank_width) {
            
            float current_plank_width = plank_width;
            if (x + plank_width > current_deck_width) {
                current_plank_width = current_deck_width - x;
            }

            generateBox({x, current_deck_y, z}, {current_plank_width, plank_height, plank_segment_length}, vertices, indices);
        }
    }
}


void generateStorageShed() {
    std::cout << "    -> Membangun kotak penyimpanan dek..." << std::endl;

    float box_length = 5.0f;
    float box_width = 1.2f;  
    float box_height = 2.2f; 
    float z_pos = 10.0f;      

    float deck_edge_x = getDeckWidthAtZ(z_pos);
    float box_x_pos = deck_edge_x - box_width - 2.0f; 

    glm::vec3 base_pos = {box_x_pos, getDeckYAtZ(z_pos), z_pos};
    
    generateBox(base_pos, {box_width, box_height, box_length}, vertices, indices);

    glm::vec3 lid_pos = base_pos + glm::vec3(-0.05f, box_height, -0.05f);
    glm::vec3 lid_size = {box_width + 0.1f, 0.1f, box_length + 0.1f};
    generateBox(lid_pos, lid_size, vertices, indices);
}


void generateAnchor() {
    glm::vec3 base_pos = {getDeckWidthAtZ(22.0f) + 0.1f, getDeckYAtZ(22.0f) - 1.0f, 22.0f};
    float shank_len = 2.0f;
    generateCylinder(base_pos, shank_len, 0.1f, 8, vertices, indices);
    generateBox({base_pos.x - 0.5f, base_pos.y, base_pos.z - 0.1f}, {1.0f, 0.2f, 0.2f}, vertices, indices);
    generateRailSegment(base_pos + glm::vec3(0, shank_len, 0), base_pos + glm::vec3(0.8f, shank_len + 0.5f, 0), 0.08f, vertices, indices);
    generateRailSegment(base_pos + glm::vec3(0, shank_len, 0), base_pos + glm::vec3(-0.8f, shank_len + 0.5f, 0), 0.08f, vertices, indices);
    base_pos.x *= -1;
    generateCylinder(base_pos, shank_len, 0.1f, 8, vertices, indices);
    generateBox({base_pos.x - 0.5f, base_pos.y, base_pos.z - 0.1f}, {1.0f, 0.2f, 0.2f}, vertices, indices);
    generateRailSegment(base_pos + glm::vec3(0, shank_len, 0), base_pos + glm::vec3(0.8f, shank_len + 0.5f, 0), 0.08f, vertices, indices);
    generateRailSegment(base_pos + glm::vec3(0, shank_len, 0), base_pos + glm::vec3(-0.8f, shank_len + 0.5f, 0), 0.08f, vertices, indices);
}

void generateAftCabin() {
    float deck_y = getDeckYAtZ(-15.0f) + 0.8f;
    float deck_width = getDeckWidthAtZ(-15.0f);

    glm::vec3 cabin_pos = {-deck_width * 0.9f, deck_y, -21.0f};
    glm::vec3 cabin_size = {deck_width * 1.8f, 0.5f, 8.0f}; 
    generateBox(cabin_pos, cabin_size, vertices, indices);

    generateBox(cabin_pos + glm::vec3(0, 0.5f, 0), {cabin_size.x, 0.6f, 0.4f}, vertices, indices);
    generateBox(cabin_pos + glm::vec3(0, 0.5f, 0), {0.4f, 0.6f, cabin_size.z}, vertices, indices);

    glm::vec3 roof_pos = cabin_pos + glm::vec3(-0.2f, 3.0f, -0.2f);
    glm::vec3 roof_size = {cabin_size.x + 0.4f, 0.2f, cabin_size.z + 0.4f};
    generateBox(roof_pos, roof_size, vertices, indices);
    generateCylinder({cabin_pos.x, deck_y, cabin_pos.z + cabin_size.z}, 3.0f, 0.1f, 12, vertices, indices);
    generateCylinder({cabin_pos.x + cabin_size.x, deck_y, cabin_pos.z + cabin_size.z}, 3.0f, 0.1f, 12, vertices, indices);

    generateBox({-1.5f, deck_y, -16.0f}, {3.0f, 0.6f, 1.5f}, vertices, indices);

    glm::vec3 wheel_pos = {0.0f, deck_y + 1.2f, -13.0f};
    generateCylinder(wheel_pos, 0.1f, 0.8f, 16, vertices, indices);
    for(int i=0; i<8; ++i) {
        float angle = 2.0f * M_PI * i / 8.0f;
        glm::vec3 spoke_dir(cos(angle), sin(angle), 0);
        generateRailSegment(wheel_pos, wheel_pos + spoke_dir * 0.7f, 0.03f, vertices, indices);
    }
}

void generatePortholes() {
    float wall_thickness = 0.2f;
    float porthole_radius = 0.3f;
    int sides = 12;

    for(float z = -10.0f; z < 15.0f; z += 4.0f) {
        float y = getDeckYAtZ(z) - 2.0f;
        float x_pos = getDeckWidthAtZ(z);
        
        glm::vec3 right_center = {x_pos, y, z};
        glm::vec3 right_orientation = {1.0f, 0.0f, 0.0f};
        glm::vec3 right_start_pos = right_center - (right_orientation * (wall_thickness / 2.0f));
        generateCylinder(right_start_pos, wall_thickness, porthole_radius, sides, vertices, indices, right_orientation);

        glm::vec3 left_center = {-x_pos, y, z};
        glm::vec3 left_orientation = {-1.0f, 0.0f, 0.0f};
        glm::vec3 left_start_pos = left_center - (left_orientation * (wall_thickness / 2.0f));
        generateCylinder(left_start_pos, wall_thickness, porthole_radius, sides, vertices, indices, left_orientation);
    }
}

void generateSunLoungers() {
    glm::vec3 pos1 = {-3.0f, getDeckYAtZ(0.0f) + 0.8f, 0.0f};
    generateBox(pos1, {2.0f, 0.3f, 4.0f}, vertices, indices);
    generateBox(pos1 + glm::vec3(0, 0.3f, 0), {2.0f, 1.0f, 0.3f}, vertices, indices);

    glm::vec3 pos2 = {1.0f, getDeckYAtZ(0.0f) + 0.8f, 0.0f};
    generateBox(pos2, {2.0f, 0.3f, 4.0f}, vertices, indices);
    generateBox(pos2 + glm::vec3(0, 0.3f, 0), {2.0f, 1.0f, 0.3f}, vertices, indices);
}


void generateSailMesh(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, int divisions, 
                      std::vector<glm::vec3>& target_vertices, 
                      std::vector<unsigned int>& target_indices) {
    int base_index = target_vertices.size();
    for (int i = 0; i <= divisions; ++i) {
        float t = static_cast<float>(i) / divisions;
        glm::vec3 side1 = glm::mix(p0, p3, t);
        glm::vec3 side2 = glm::mix(p1, p2, t);
        for (int j = 0; j <= divisions; ++j) {
            float s = static_cast<float>(j) / divisions;
            glm::vec3 point = glm::mix(side1, side2, s);
            
            if (glm::length(p1-p0) > 0.1f && glm::length(p3-p0) > 0.1f) {
                glm::vec3 sail_normal = glm::normalize(glm::cross(p1 - p0, p3 - p0));
                float billow = sin(t * M_PI) * sin(s * M_PI) * 2.0f;
                point += sail_normal * billow;
            }
            target_vertices.push_back(point);
        }
    }
    for (int i = 0; i < divisions; ++i) {
        for (int j = 0; j < divisions; ++j) {
            int row1 = i * (divisions + 1);
            int row2 = (i + 1) * (divisions + 1);
            target_indices.push_back(base_index + row1 + j);
            target_indices.push_back(base_index + row2 + j + 1);
            target_indices.push_back(base_index + row1 + j + 1);
            
            target_indices.push_back(base_index + row1 + j);
            target_indices.push_back(base_index + row2 + j);
            target_indices.push_back(base_index + row2 + j + 1);
        }
    }
}

float interpolate(float val, const std::vector<float>& x, const std::vector<float>& y) {
    if (val <= x.front()) return y.front();
    if (val >= x.back()) return y.back();
    auto it = std::upper_bound(x.begin(), x.end(), val);
    int idx = std::distance(x.begin(), it) - 1;
    float t = (val - x[idx]) / (x[idx+1] - x[idx]);
    return glm::mix(y[idx], y[idx+1], t);
}

float getDeckYAtZ(float z) {
    std::vector<float> y_points;
    for(const auto& section : blueprint_sections) {
        y_points.push_back(section.getPoint(1.0f).y);
    }
    return interpolate(z, blueprint_z_positions, y_points);
}

float getDeckWidthAtZ(float z) {
    std::vector<float> x_points;
    for(const auto& section : blueprint_sections) {
        x_points.push_back(section.getPoint(1.0f).x);
    }
    return interpolate(z, blueprint_z_positions, x_points);
}


// FUNGSI OPENGL, KONTROL & MAIN

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 1280.0/960.0, 1.0, 200.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    float rad_x = camera_angle_x * M_PI / 180.0f;
    float rad_y = camera_angle_y * M_PI / 180.0f;
    float eye_offset_x = camera_distance * cos(rad_x) * sin(rad_y);
    float eye_offset_y = camera_distance * sin(rad_x);
    float eye_offset_z = camera_distance * cos(rad_x) * cos(rad_y);
    glm::vec3 eye_pos = camera_target + glm::vec3(eye_offset_x, eye_offset_y, eye_offset_z);
    gluLookAt(eye_pos.x, eye_pos.y, eye_pos.z, camera_target.x, camera_target.y, camera_target.z, 0.0, 1.0, 0.0);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_NORMALIZE);

    GLfloat light0_pos[] = {-50.0f, 50.0f, 50.0f, 1.0f};
    GLfloat light0_diffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);

    GLfloat light1_pos[] = {50.0f, 20.0f, -50.0f, 1.0f};
    GLfloat light1_diffuse[] = {0.4f, 0.4f, 0.4f, 1.0f};
    glLightfv(GL_LIGHT1, GL_POSITION, light1_pos);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, wood_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, wood_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, wood_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, wood_shininess);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, vertices.data());
    glNormalPointer(GL_FLOAT, 0, normals.data());
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, indices.data());
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, sail_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, sail_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, sail_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, sail_shininess);

    glVertexPointer(3, GL_FLOAT, 0, sail_vertices.data());
    glNormalPointer(GL_FLOAT, 0, sail_normals.data());
    glDrawElements(GL_TRIANGLES, sail_indices.size(), GL_UNSIGNED_INT, sail_indices.data());

    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);

    glutSwapBuffers();
}

void mouseButton(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        is_mouse_dragging = (state == GLUT_DOWN);
        last_mouse_x = x;
        last_mouse_y = y;
    }
}

void mouseMove(int x, int y) {
    if (is_mouse_dragging) {
        camera_angle_y += (x - last_mouse_x) * 0.4f;
        camera_angle_x += (y - last_mouse_y) * 0.4f;
        camera_angle_x = glm::clamp(camera_angle_x, -89.0f, 89.0f);
        last_mouse_x = x;
        last_mouse_y = y;
        glutPostRedisplay();
    }
}

void mouseWheel(int button, int dir, int x, int y) {
    camera_distance -= dir * 2.0f;
    camera_distance = glm::clamp(camera_distance, 10.0f, 200.0f);
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
    float pan_speed = 1.0f;
    if (key == 'w') camera_target.z -= pan_speed;
    if (key == 's') camera_target.z += pan_speed;
    if (key == 'a') camera_target.x -= pan_speed;
    if (key == 'd') camera_target.x += pan_speed;
    if (key == 'q') camera_target.y += pan_speed;
    if (key == 'e') camera_target.y -= pan_speed;
    if (key == 'r') {
        camera_distance = 80.0f;
        camera_angle_x = 30.0f;
        camera_angle_y = -60.0f;
        camera_target = {0.0f, 15.0f, 0.0f};
    }
    glutPostRedisplay();
}

void exportToObj(const char* filename) {
    std::ofstream objFile(filename);
    if (!objFile.is_open()) {
        std::cerr << "Error: Tidak bisa membuka file untuk ekspor!" << std::endl;
        return;
    }

    objFile << "# Model Kapal Phinisi Modern v2.8 - The Grand Finale" << std::endl;
    objFile << "o Phinisi_Schooner" << std::endl;

    for (const auto& v : vertices) objFile << "v " << v.x << " " << v.y << " " << v.z << std::endl;
    for (const auto& n : normals) objFile << "vn " << n.x << " " << n.y << " " << n.z << std::endl;
    
    objFile << "usemtl Wood" << std::endl;
    for (size_t i = 0; i < indices.size(); i += 3) {
        objFile << "f " << indices[i]+1 << "//" << indices[i]+1 << " " << indices[i+1]+1 << "//" << indices[i+1]+1 << " " << indices[i+2]+1 << "//" << indices[i+2]+1 << std::endl;
    }

    int vertex_offset = vertices.size();
    for (const auto& v : sail_vertices) objFile << "v " << v.x << " " << v.y << " " << v.z << std::endl;
    for (const auto& n : sail_normals) objFile << "vn " << n.x << " " << n.y << " " << n.z << std::endl;

    objFile << "usemtl Sail" << std::endl;
    for (size_t i = 0; i < sail_indices.size(); i += 3) {
        int i1 = vertex_offset + sail_indices[i] + 1;
        int i2 = vertex_offset + sail_indices[i+1] + 1;
        int i3 = vertex_offset + sail_indices[i+2] + 1;
        objFile << "f " << i1 << "//" << i1 << " " << i2 << "//" << i2 << " " << i3 << "//" << i3 << std::endl;
    }

    objFile.close();
    std::cout << "Model berhasil diekspor ke " << filename << std::endl;
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    srand(time(NULL));
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1280, 960);
    glutCreateWindow("Kapal Pinisi - kelompok 2");
    
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    GLfloat ambient_light[] = {0.4f, 0.4f, 0.4f, 1.0f}; 
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);
    
    std::cout << "Memulai rekayasa kapal Phinisi modern v2.8..." << std::endl;
    generateHullAndDeck();
    std::cout << "Lambung selesai." << std::endl;
    generateAftCabin();
    std::cout << "Kabin mewah belakang terpasang." << std::endl;
    generateDeckPlanking();
    std::cout << "Papan dek terpasang." << std::endl;
    std::cout << "Platform bowsprit kompleks sedang dibangun..." << std::endl;
    glm::vec3 bowsprit_tip = generateBowspritPlatform(); 
    std::cout << "Platform bowsprit selesai." << std::endl;
    generateSideNetting(true, vertices, indices); 
    generateSideNetting(false, vertices, indices);
    generateSails();
    std::cout << "Layar terkembang." << std::endl;
    std::cout << "Tiang dan spar sedang dibuat..." << std::endl;
    MastGeometry foremast_geo = generateSingleMast({0.0f, getDeckYAtZ(12.0f), 12.0f}, 45.0f, true);
    MastGeometry mainmast_geo = generateSingleMast({0.0f, getDeckYAtZ(-12.0f), -12.0f}, 40.0f, false);
    std::cout << "Tiang dan spar terpasang." << std::endl;
    std::cout << "Tali-temali sedang dipasang secara presisi..." << std::endl;
    generateRigging_Phinisi(foremast_geo, mainmast_geo, bowsprit_tip);
    std::cout << "Tali-temali presisi selesai." << std::endl;
    generateMastLadder(foremast_geo.base, 45.0f, 0.45f, {-1.0f, 0.0f, 0.0f}, vertices, indices);
    
    generateMastLadder(mainmast_geo.base, 40.0f, 0.45f, {0.0f, 0.0f, 1.0f}, vertices, indices);
    generateRudder();
    std::cout << "Kemudi terpasang." << std::endl;
    generateDeckRailings();
    std::cout << "Pagar dek dan anjungan depan selesai." << std::endl;
    generateAnchor();
    std::cout << "Jangkar terpasang." << std::endl;
    generatePortholes();
    std::cout << "Jendela lambung terpasang." << std::endl;
    generateDeckhouse();
    generateRooftopDining();
    generateBarrelCluster({0.0f, getDeckYAtZ(-8.0f), -8.0f});
    generateBarrelCluster({-2.5f, getDeckYAtZ(10.0f), 10.0f});
    generateStorageShed();
    placeCoiledRopes();
    placeLifebuoys();
    generateSunLoungers();
    std::cout << "Perabotan mewah ditambahkan." << std::endl;

    recalculateNormals(vertices, indices, normals);
    recalculateNormals(sail_vertices, sail_indices, sail_normals);
    std::cout << "Pencahayaan dikalkulasi ulang. Model siap." << std::endl;
    
    exportToObj("pinisi3D.obj");

    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMove);
    glutMouseWheelFunc(mouseWheel);
    
    glutMainLoop();
    
    return 0;
}