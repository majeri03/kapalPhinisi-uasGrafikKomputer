#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iterator>

// DEFINISI WAJIB UNTUK MENGGUNAKAN FUNGSI EKSPERIMENTAL GLM
#define GLM_ENABLE_EXPERIMENTAL

#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ===================================================================
// DEKLARASI FUNGSI (FORWARD DECLARATION)
// ===================================================================
void generateSailMesh(glm::vec3 p0, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, int divisions, std::vector<glm::vec3>& target_vertices, std::vector<unsigned int>& target_indices);
float getDeckYAtZ(float z);
float getDeckWidthAtZ(float z);

// ===================================================================
// VARIABEL GLOBAL & KONFIGURASI
// ===================================================================

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

// REVISI: Warna layar diubah menjadi biru tua
GLfloat sail_ambient[] = {0.05f, 0.1f, 0.2f, 1.0f};
GLfloat sail_diffuse[] = {0.1f, 0.2f, 0.4f, 1.0f};
GLfloat sail_specular[] = {0.1f, 0.1f, 0.1f, 1.0f};
GLfloat sail_shininess[] = {30.0f};

// ===================================================================
// STRUKTUR & BLUEPRINT DASAR
// ===================================================================

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


// ===================================================================
// FUNGSI-FUNGSI GENERATOR GEOMETRI
// ===================================================================

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

void generateCylinder(glm::vec3 pos, float height, float radius, int sides, std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref) {
    int base_index = vert_ref.size();
    vert_ref.push_back(pos);
    vert_ref.push_back(pos + glm::vec3(0, height, 0));

    for (int i = 0; i <= sides; ++i) {
        float angle = 2.0f * M_PI * i / sides;
        glm::vec3 offset(cos(angle) * radius, 0, sin(angle) * radius);
        vert_ref.push_back(pos + offset);
        vert_ref.push_back(pos + glm::vec3(0, height, 0) + offset);
    }

    for (int i = 0; i < sides; ++i) {
        int i0 = base_index + 2 + i * 2;
        int i1 = base_index + 2 + (i + 1) * 2;
        ind_ref.push_back(i0); ind_ref.push_back(i1); ind_ref.push_back(i0 + 1);
        ind_ref.push_back(i1); ind_ref.push_back(i1 + 1); ind_ref.push_back(i0 + 1);
    }
    
    for (int i = 0; i < sides; ++i) {
        int i0 = base_index + 2 + i * 2;
        int i1 = base_index + 2 + (i + 1) * 2;
        ind_ref.push_back(base_index); ind_ref.push_back(i0); ind_ref.push_back(i1);
        ind_ref.push_back(base_index + 1); ind_ref.push_back(i1 + 1); ind_ref.push_back(i0 + 1);
    }
}

void generateRailSegment(glm::vec3 p1, glm::vec3 p2, float radius, std::vector<glm::vec3>& vert_ref, std::vector<unsigned int>& ind_ref) {
    glm::vec3 dir = p2 - p1;
    if (glm::length(dir) < 0.001f) return;
    glm::quat rot = glm::rotation(glm::vec3(0, 1, 0), glm::normalize(dir));

    int base_index = vert_ref.size();
    int sides = 6;

    for (int i = 0; i <= sides; ++i) {
        float angle = 2.0f * M_PI * i / sides;
        glm::vec3 offset(cos(angle) * radius, 0, sin(angle) * radius);
        offset = rot * offset;
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

void generateMastsAndSpars() {
    // REVISI: Tiang-tiang ditinggikan secara signifikan
    glm::vec3 foremast_pos = {0.0f, getDeckYAtZ(12.0f), 12.0f};
    float foremast_height = 45.0f; // Ditinggikan
    glm::vec3 mainmast_pos = {0.0f, getDeckYAtZ(-12.0f), -12.0f};
    float mainmast_height = 40.0f; // Ditinggikan
    
    generateCylinder(foremast_pos, foremast_height, 0.45f, 16, vertices, indices);
    generateCylinder(mainmast_pos, mainmast_height, 0.4f, 16, vertices, indices);

    // Spreaders disesuaikan
    generateBox({-3.0f, foremast_pos.y + 22.0f, foremast_pos.z - 0.2f}, {6.0f, 0.3f, 0.4f}, vertices, indices);
    generateBox({-2.5f, mainmast_pos.y + 20.0f, mainmast_pos.z - 0.2f}, {5.0f, 0.3f, 0.4f}, vertices, indices);

    glm::vec3 bowsprit_base = {0.0f, getDeckYAtZ(25.0f) + 0.8f, 25.0f};
    glm::vec3 bowsprit_tip = {0.0f, bowsprit_base.y + 2.0f, 40.0f};
    generateRailSegment(bowsprit_base, bowsprit_tip, 0.3f, vertices, indices);
    generateRailSegment(bowsprit_tip, {0.0f, getDeckYAtZ(25.0f) - 2.0f, 26.0f}, 0.05f, vertices, indices);
    generateRailSegment({0.0f, getDeckYAtZ(25.0f) - 2.0f, 26.0f}, {2.0f, getDeckYAtZ(24.0f), 24.0f}, 0.05f, vertices, indices);
    generateRailSegment({0.0f, getDeckYAtZ(25.0f) - 2.0f, 26.0f}, {-2.0f, getDeckYAtZ(24.0f), 24.0f}, 0.05f, vertices, indices);
    
    // Lampu Navigasi
    generateCylinder(foremast_pos + glm::vec3(0, foremast_height, 0), 0.3f, 0.15f, 8, vertices, indices);
    generateCylinder(mainmast_pos + glm::vec3(0, mainmast_height, 0), 0.3f, 0.15f, 8, vertices, indices);
}

void generateSails() {
    glm::vec3 mainmast_pos = {0.0f, getDeckYAtZ(-12.0f), -12.0f};
    glm::vec3 foremast_pos = {0.0f, getDeckYAtZ(12.0f), 12.0f};
    glm::vec3 bowsprit_tip = {0.0f, getDeckYAtZ(25.0f) + 2.8f, 40.0f};

    // REVISI: Ketinggian puncak tiang disesuaikan
    glm::vec3 mainmast_top = mainmast_pos + glm::vec3(0, 40.0f, 0);
    glm::vec3 foremast_top = foremast_pos + glm::vec3(0, 45.0f, 0);

    // REVISI: Layar Utama Belakang (Segitiga), posisi boom dinaikkan
    glm::vec3 main_p0 = mainmast_top;
    glm::vec3 main_p1 = mainmast_pos + glm::vec3(0, 4.0f, 0); // Dinaikkan agar tidak menimpa kabin
    glm::vec3 main_p2 = {0.0f, main_p1.y, -23.0f};
    generateRailSegment(main_p1, main_p2, 0.25f, vertices, indices); // Boom
    generateSailMesh(main_p0, main_p1, main_p2, main_p2, 15, sail_vertices, sail_indices);

    // REVISI: Satu layar segitiga di belakang tiang depan (Foresail)
    glm::vec3 fore_p0 = foremast_top;
    glm::vec3 fore_p1 = foremast_pos + glm::vec3(0, 1.5f, 0);
    glm::vec3 fore_p2 = {0.0f, fore_p1.y, -2.0f};
    generateRailSegment(fore_p1, fore_p2, 0.25f, vertices, indices); // Boom
    generateSailMesh(fore_p0, fore_p1, fore_p2, fore_p2, 15, sail_vertices, sail_indices);

    // REVISI: Dua Layar Jib di Depan Tiang Depan
    glm::vec3 jib_clew = {0.0f, getDeckYAtZ(20.0f) + 1.0f, 20.0f};
    generateSailMesh(foremast_top, jib_clew, bowsprit_tip, bowsprit_tip, 12, sail_vertices, sail_indices);
    generateSailMesh({foremast_top.x, foremast_top.y - 15.0f, foremast_pos.z}, jib_clew, {bowsprit_tip.x, bowsprit_tip.y-1.0f, bowsprit_tip.z - 12.0f}, {bowsprit_tip.x, bowsprit_tip.y-1.0f, bowsprit_tip.z - 12.0f}, 12, sail_vertices, sail_indices);
}

void generateRiggingAndDetails() {
    float rope_radius = 0.05f;
    float thin_rope_radius = 0.03f;

    glm::vec3 mainmast_pos = {0.0f, getDeckYAtZ(-12.0f), -12.0f};
    glm::vec3 foremast_pos = {0.0f, getDeckYAtZ(12.0f), 12.0f};
    glm::vec3 bowsprit_tip = {0.0f, getDeckYAtZ(25.0f) + 2.8f, 40.0f};
    glm::vec3 mainmast_top = mainmast_pos + glm::vec3(0, 40.0f, 0);
    glm::vec3 foremast_top = foremast_pos + glm::vec3(0, 45.0f, 0);

    // Tali Utama
    generateRailSegment(foremast_top, bowsprit_tip, rope_radius, vertices, indices);
    generateRailSegment(mainmast_top, foremast_top, rope_radius, vertices, indices);
    generateRailSegment(mainmast_top, {0.0f, getDeckYAtZ(-25.0f), -25.0f}, rope_radius, vertices, indices);

    // Tali Kompleks dengan Spreader
    auto generateComplexShrouds = [&](glm::vec3 mast_pos, float mast_height, float deck_y, float z_pos, float width, float spreader_y, float spreader_width) {
        glm::vec3 spreader_left = {-spreader_width/2, spreader_y, z_pos};
        glm::vec3 spreader_right = {spreader_width/2, spreader_y, z_pos};
        glm::vec3 deck_left = {-width/2, deck_y, z_pos};
        glm::vec3 deck_right = {width/2, deck_y, z_pos};
        glm::vec3 mast_top = mast_pos + glm::vec3(0, mast_height, 0);

        generateRailSegment(mast_top, spreader_left, thin_rope_radius, vertices, indices);
        generateRailSegment(mast_top, spreader_right, thin_rope_radius, vertices, indices);
        generateRailSegment(spreader_left, deck_left, thin_rope_radius, vertices, indices);
        generateRailSegment(spreader_right, deck_right, thin_rope_radius, vertices, indices);
        generateRailSegment(mast_pos + glm::vec3(0, spreader_y, 0), deck_left, thin_rope_radius, vertices, indices);
        generateRailSegment(mast_pos + glm::vec3(0, spreader_y, 0), deck_right, thin_rope_radius, vertices, indices);
    };

    generateComplexShrouds(foremast_pos, 45.0f, getDeckYAtZ(12.0f) + 0.8f, 12.0f, 8.0f, foremast_pos.y + 22.0f, 6.0f);
    generateComplexShrouds(mainmast_pos, 40.0f, getDeckYAtZ(-12.0f) + 0.8f, -12.0f, 7.5f, mainmast_pos.y + 20.0f, 5.0f);
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
    // Bow Pulpit (Pagar Depan)
    generateRailSegment(right_posts.back(), left_posts.back(), post_radius, vertices, indices);
}

void generateDeckPlanking() {
    float plank_width = 0.2f;
    float plank_height = 0.05f;
    for(float x = -getDeckWidthAtZ(0.0f) + plank_width; x < getDeckWidthAtZ(0.0f); x += plank_width) {
        generateBox({x, getDeckYAtZ(-20.0f) - plank_height, -20.0f}, {plank_width, plank_height, 45.0f}, vertices, indices);
    }
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
    glm::vec3 cabin_size = {deck_width * 1.8f, 0.5f, 8.0f}; // Base for sofa
    generateBox(cabin_pos, cabin_size, vertices, indices);

    // Sandaran Sofa
    generateBox(cabin_pos + glm::vec3(0, 0.5f, 0), {cabin_size.x, 0.6f, 0.4f}, vertices, indices);
    generateBox(cabin_pos + glm::vec3(0, 0.5f, 0), {0.4f, 0.6f, cabin_size.z}, vertices, indices);

    // Atap & Pilar
    glm::vec3 roof_pos = cabin_pos + glm::vec3(-0.2f, 3.0f, -0.2f);
    glm::vec3 roof_size = {cabin_size.x + 0.4f, 0.2f, cabin_size.z + 0.4f};
    generateBox(roof_pos, roof_size, vertices, indices);
    generateCylinder({cabin_pos.x, deck_y, cabin_pos.z + cabin_size.z}, 3.0f, 0.1f, 12, vertices, indices);
    generateCylinder({cabin_pos.x + cabin_size.x, deck_y, cabin_pos.z + cabin_size.z}, 3.0f, 0.1f, 12, vertices, indices);

    // Meja
    generateBox({-1.5f, deck_y, -16.0f}, {3.0f, 0.6f, 1.5f}, vertices, indices);

    // Roda Kemudi
    glm::vec3 wheel_pos = {0.0f, deck_y + 1.2f, -13.0f};
    generateCylinder(wheel_pos, 0.1f, 0.8f, 16, vertices, indices);
    for(int i=0; i<8; ++i) {
        float angle = 2.0f * M_PI * i / 8.0f;
        glm::vec3 spoke_dir(cos(angle), sin(angle), 0);
        generateRailSegment(wheel_pos, wheel_pos + spoke_dir * 0.7f, 0.03f, vertices, indices);
    }
}

void generatePortholes() {
    for(float z = -10.0f; z < 15.0f; z += 4.0f) {
        float y = getDeckYAtZ(z) - 2.0f;
        float x = getDeckWidthAtZ(z) + 0.05f;
        generateCylinder({x, y, z}, 0.1f, 0.3f, 12, vertices, indices);
        generateCylinder({-x, y, z}, 0.1f, 0.3f, 12, vertices, indices);
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


// ===================================================================
// FUNGSI OPENGL, KONTROL & MAIN
// ===================================================================

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
    glEnable(GL_LIGHT1); // Aktifkan cahaya kedua
    glEnable(GL_NORMALIZE);

    // Cahaya Utama (Key Light)
    GLfloat light0_pos[] = {-50.0f, 50.0f, 50.0f, 1.0f};
    GLfloat light0_diffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);

    // Cahaya Pengisi (Fill Light) untuk mengurangi bayangan
    GLfloat light1_pos[] = {50.0f, 20.0f, -50.0f, 1.0f};
    GLfloat light1_diffuse[] = {0.4f, 0.4f, 0.4f, 1.0f}; // Lebih redup
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
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1280, 960);
    glutCreateWindow("Kapal Phinisi Modern - v2.8 The Grand Finale");
    
    // REVISI: Latar belakang diubah menjadi putih
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    GLfloat ambient_light[] = {0.4f, 0.4f, 0.4f, 1.0f}; // Cahaya ambien sedikit lebih terang
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light);
    
    std::cout << "Memulai rekayasa kapal Phinisi modern v2.8..." << std::endl;
    generateHullAndDeck();
    std::cout << "Lambung selesai." << std::endl;
    generateAftCabin();
    std::cout << "Kabin mewah belakang terpasang." << std::endl;
    generateDeckPlanking();
    std::cout << "Papan dek terpasang." << std::endl;
    generateMastsAndSpars();
    std::cout << "Tiang dan spar terpasang." << std::endl;
    generateSails();
    std::cout << "Layar terkembang." << std::endl;
    generateRiggingAndDetails();
    std::cout << "Tali-temali kompleks selesai." << std::endl;
    generateRudder();
    std::cout << "Kemudi terpasang." << std::endl;
    generateDeckRailings();
    std::cout << "Pagar dek dan anjungan depan selesai." << std::endl;
    generateAnchor();
    std::cout << "Jangkar terpasang." << std::endl;
    generatePortholes();
    std::cout << "Jendela lambung terpasang." << std::endl;
    generateSunLoungers();
    std::cout << "Perabotan mewah ditambahkan." << std::endl;

    recalculateNormals(vertices, indices, normals);
    recalculateNormals(sail_vertices, sail_indices, sail_normals);
    std::cout << "Pencahayaan dikalkulasi ulang. Model siap." << std::endl;
    
    exportToObj("phinisi_modern_v2.8.obj");

    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMove);
    glutMouseWheelFunc(mouseWheel);
    
    glutMainLoop();
    
    return 0;
}
