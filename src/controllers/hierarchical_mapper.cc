// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "controllers/hierarchical_mapper.h"

#include "base/scene_clustering.h"
#include "util/misc.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

namespace colmap {
namespace {

// cluster为引用，reconstruction_managers为指针
void MergeClusters(
    const SceneClustering::Cluster& cluster,
    std::unordered_map<const SceneClustering::Cluster*, ReconstructionManager>*
        reconstruction_managers) {
  std::cout << " 进入Merge " << std::endl;
  std::cout << " &cluster = " << &cluster << std::endl;
  std::cout << " cluster.child_clusters.size() = " << cluster.child_clusters.size() << std::endl;

  for (size_t i = 0; i < cluster.child_clusters.size(); i++)
  {
    std::cout << "cluster.child_cluster[ " << i << " ] = " << &cluster.child_clusters[i] << std::endl;
  }
  

  // FIXME reconstruction_managers = cluster 的数目
  std::cout << "reconstruction_managers = " << reconstruction_managers << std::endl;  

  std::cout << "reconstruction_managers.size() = " << reconstruction_managers->size() << std::endl;

  // FIXME 每个cluster对应一个reconstruction
  for (const auto& recon_man : *reconstruction_managers)
  {
    std::cout << " recon_man.first = " << recon_man.first <<  " recon_man.first->image_ids.size() = " << 
    recon_man.first->image_ids.size() << " , recon_man.second = " << recon_man.second.Size() << std::endl;
  }
  

  std::vector<Reconstruction*> reconstructions;
  for (const auto& child_cluster : cluster.child_clusters) {
    std::cout << " Merge child_cluster " << std::endl;
    if (!child_cluster.child_clusters.empty()) {  //孩子节点的孩子节点不为空，执行merge
      std::cout << " 孩子节点的孩子节点不为空，执行merge  " << std::endl;
      MergeClusters(child_cluster, reconstruction_managers);
    }

    std::cout << "child_cluster.image_ids.size() = " << child_cluster.image_ids.size() << std::endl;
    std::cout << "child_cluster的地址为 ：" << &child_cluster << std::endl;

    std::cout << "reconstruction_managers->size = " << reconstruction_managers->size() << std::endl;
    // std::cout << "reconstruction_manager->first = " << reconstruction_managers->first() << std::endl;

    // FIXME at返回value值，即返回的数据类型是ReconstructionManager
    auto& reconstruction_manager = reconstruction_managers->at(&child_cluster);
    
    std::cout << " reconstruction_manager.Size() = " << reconstruction_manager.Size() << std::endl;

    for (size_t i = 0; i < reconstruction_manager.Size(); ++i) {  //将孩子节点重建的结果放入容器reconstructions中
      std::cout << "将重建结果放入reconstructions " << std::endl;
      // FIXME 将重建结果 reconstruction 放入 reconstructions中
      reconstructions.push_back(&reconstruction_manager.Get(i));
    }
  }

  std::cout << " reconstructions.size() =  " << reconstructions.size() << std::endl;

  // Try to merge all child cluster reconstruction.
  while (reconstructions.size() > 1) {
    std::cout << " 开始meige child clusters " << std::endl;
    bool merge_success = false;
    for (size_t i = 0; i < reconstructions.size(); ++i) {
      for (size_t j = 0; j < i; ++j) {
        const double kMaxReprojError = 8.0;
        if (reconstructions[i]->Merge(*reconstructions[j], kMaxReprojError)) {
          reconstructions.erase(reconstructions.begin() + j);
          merge_success = true;
          break;
        }
      }

      if (merge_success) {
        std::cout << "merge success" << std::endl;
        break;
      }
    }

    if (!merge_success) {
      std::cout << "merge failed" << std::endl;
      break;
    }
  }

  //FIXME unordered_map 中 [] 相当于at，在unordered_map中寻找是不是存在cluster
  // * 取出指针变量reconstruction_managers中的值
  // Insert a new reconstruction manager for merged cluster.
  auto& reconstruction_manager = (*reconstruction_managers)[&cluster];
  for (const auto& reconstruction : reconstructions) {
    std::cout << " reconstruction_manager.Size = " << reconstruction_manager.Size() << std::endl;
    std::cout << " Insert a new reconstruction manager " << std::endl;
    reconstruction_manager.Add();
    reconstruction_manager.Get(reconstruction_manager.Size() - 1) =
        *reconstruction;
  }

  // Delete all merged child cluster reconstruction managers.
  for (const auto& child_cluster : cluster.child_clusters) {
    std::cout << "delete merged" << std::endl;
    reconstruction_managers->erase(&child_cluster);
  }
}

}  // namespace

bool HierarchicalMapperController::Options::Check() const {
  CHECK_OPTION_GT(init_num_trials, -1);
  CHECK_OPTION_GE(num_workers, -1);
  return true;
}

HierarchicalMapperController::HierarchicalMapperController(
    const Options& options, const SceneClustering::Options& clustering_options,
    const IncrementalMapperOptions& mapper_options,
    ReconstructionManager* reconstruction_manager)
    : options_(options),
      clustering_options_(clustering_options),
      mapper_options_(mapper_options),
      reconstruction_manager_(reconstruction_manager) {
  CHECK(options_.Check());
  CHECK(clustering_options_.Check());
  CHECK(mapper_options_.Check());
  CHECK_EQ(clustering_options_.branching, 2);
}

void HierarchicalMapperController::Run() {
  PrintHeading1("Partitioning the scene");

  //////////////////////////////////////////////////////////////////////////////
  // Cluster scene
  //////////////////////////////////////////////////////////////////////////////

  SceneClustering scene_clustering(clustering_options_);

  std::unordered_map<image_t, std::string> image_id_to_name;

  {
    Database database(options_.database_path);

    std::cout << "Reading images..." << std::endl;
    const auto images = database.ReadAllImages();
    for (const auto& image : images) {
      image_id_to_name.emplace(image.ImageId(), image.Name());
    }

    std::cout << "Reading scene graph..." << std::endl;
    std::vector<std::pair<image_t, image_t>> image_pairs;
    std::vector<int> num_inliers;
    database.ReadTwoViewGeometryNumInliers(&image_pairs, &num_inliers);

    std::cout << "Partitioning scene graph..." << std::endl;
    // scene_clustering.Partition(image_pairs, num_inliers);
  }


  //FIXME 源代码 分好组的image_id  
  // auto leaf_clusters = scene_clustering.GetLeafClusters();
  
  //FIXME 读入自己分组的结果进行重建
  // 按行读入
  std::vector<const colmap::SceneClustering::Cluster *> leaf_clusters;
  colmap::SceneClustering::Cluster * const c1 = new colmap::SceneClustering::Cluster();
  std::ifstream c1txt("/data6/wcm/data/Building/NCut/4/cluster_id/overlap/c0.txt");
  std::string w;
  while (getline(c1txt, w))
  {
    std::stringstream c1txt(w);
    std::string tmp;
    while (getline(c1txt,tmp))
    {
      int x = std::stoi(tmp);
      colmap::image_t t = x;
      c1->image_ids.push_back(t);
    }  
  }
  std::cout << "c1->image_ids.size = " << c1->image_ids.size() << std::endl;

  colmap::SceneClustering::Cluster * const c0 = new colmap::SceneClustering::Cluster();
  std::ifstream c0txt("/data6/wcm/data/Building/NCut/4/cluster_id/overlap/c1.txt");
  std::string w0;
  while (getline(c0txt, w0))
  {
    std::stringstream c0txt(w0);
    std::string tmp0;
    while (getline(c0txt,tmp0))
    {
      int x = std::stoi(tmp0);
      colmap::image_t t = x;
      c0->image_ids.push_back(t);
    }  
  }
 std::cout << "c0->image_ids.size = " << c0->image_ids.size() << std::endl;

  colmap::SceneClustering::Cluster * const c2 = new colmap::SceneClustering::Cluster();
  std::ifstream c2txt("/data6/wcm/data/Building/NCut/4/cluster_id/overlap/c2.txt");
  std::string w2;
  while (getline(c2txt, w2))
  {
    std::stringstream c2txt(w2);
    std::string tmp2;
    while (getline(c2txt,tmp2))
    {
      int x = std::stoi(tmp2);
      colmap::image_t t = x;
      c2->image_ids.push_back(t);
    }  
  }
 std::cout << "c2->image_ids.size = " << c2->image_ids.size() << std::endl;

  colmap::SceneClustering::Cluster * const c3 = new colmap::SceneClustering::Cluster();
  std::ifstream c3txt("/data6/wcm/data/Building/NCut/4/cluster_id/overlap/c3.txt");
  std::string w3;
  while (getline(c3txt, w3))
  {
    std::stringstream c3txt(w3);
    std::string tmp3;
    while (getline(c3txt,tmp3))
    {
      int x = std::stoi(tmp3);
      colmap::image_t t = x;
      c3->image_ids.push_back(t);
    }  
  }
 std::cout << "c3->image_ids.size = " << c3->image_ids.size() << std::endl;

  colmap::SceneClustering::Cluster * const c4 = new colmap::SceneClustering::Cluster();
  std::ifstream c4txt("/data6/wcm/data/Building/NCut/4/cluster_id/overlap/c4.txt");
  std::string w4;
  while (getline(c4txt, w4))
  {
    std::stringstream c4txt(w4);
    std::string tmp4;
    while (getline(c4txt,tmp4))
    {
      int x = std::stoi(tmp4);
      colmap::image_t t = x;
      c4->image_ids.push_back(t);
    }  
  }
 std::cout << "c4->image_ids.size = " << c4->image_ids.size() << std::endl;

 leaf_clusters.push_back(c0);
 leaf_clusters.push_back(c1);
 leaf_clusters.push_back(c2);
 leaf_clusters.push_back(c3);
 leaf_clusters.push_back(c4);
 
//FIXME 
 scene_clustering.GetLeafClusters() = leaf_clusters;

  // FIXME 源代码
  size_t total_num_images = 0;
  for (size_t i = 0; i < leaf_clusters.size(); ++i) {
    total_num_images += leaf_clusters[i]->image_ids.size();
    std::cout << StringPrintf("  Cluster %d with %d images", i + 1,
                              leaf_clusters[i]->image_ids.size())
              << std::endl;
    // std::cout << "i = " << i << ", image_id = " << std::endl;
    // for(size_t j = 0;j < leaf_clusters[i]->image_ids.size();j ++){
    //   std::cout << leaf_clusters[i]->image_ids[j] << std::endl;
    // }

  }

 

  std::cout << StringPrintf("Clusters have %d images", total_num_images)
            << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  // Reconstruct clusters
  //////////////////////////////////////////////////////////////////////////////

  PrintHeading1("Reconstructing clusters");

  // // Determine the number of workers and threads per worker.
  const int kMaxNumThreads = -1;
  const int num_eff_threads = GetEffectiveNumThreads(kMaxNumThreads);
  const int kDefaultNumWorkers = 8;
  const int num_eff_workers =
      options_.num_workers < 1
          ? std::min(static_cast<int>(leaf_clusters.size()),
                     std::min(kDefaultNumWorkers, num_eff_threads))
          : options_.num_workers;
  const int num_threads_per_worker =
      std::max(1, num_eff_threads / num_eff_workers);

  // Function to reconstruct one cluster using incremental mapping.
  auto ReconstructCluster = [&, this](
                                const SceneClustering::Cluster& cluster,
                                ReconstructionManager* reconstruction_manager) {
    if (cluster.image_ids.empty()) {
      return;
    }

    //FIXME 检查reconstruction_manager的大小
    std::cout << "reconstruction_manager.size() = " << reconstruction_manager->Size() << std::endl;

    IncrementalMapperOptions custom_options = mapper_options_;
    custom_options.max_model_overlap = 3;
    custom_options.init_num_trials = options_.init_num_trials;
    custom_options.num_threads = num_threads_per_worker;

    for (const auto image_id : cluster.image_ids) {
      custom_options.image_names.insert(image_id_to_name.at(image_id));
    } 

    IncrementalMapperController mapper(&custom_options, options_.image_path,
                                       options_.database_path,
                                       reconstruction_manager);
    mapper.Start();
    mapper.Wait();
  };

  // Start reconstructing the bigger clusters first for resource usage.
  std::sort(leaf_clusters.begin(), leaf_clusters.end(),
            [](const SceneClustering::Cluster* cluster1,
               const SceneClustering::Cluster* cluster2) {
              return cluster1->image_ids.size() > cluster2->image_ids.size();
            });

  // Start the reconstruction workers.

  std::unordered_map<const SceneClustering::Cluster*, ReconstructionManager>
      reconstruction_managers;
  // reconstruction_managers.reserve(scene_clustering.GetLeafClusters().size()); //reserve开辟空间
  reconstruction_managers.reserve(leaf_clusters.size()); //reserve开辟空间

  std::cout << "scene_clustering.GetLeafClusters().size() = " << scene_clustering.GetLeafClusters().size() << std::endl;
  std::cout << "reconstruction_managers.size = " << reconstruction_managers.size() << std::endl;

  ThreadPool thread_pool(num_eff_workers);

  // for (const auto& cluster : scene_clustering.GetLeafClusters())  
  // for (const auto& cluster : leaf_clusters)
  // {
  //   //FIXME 开始重建工作
  //   auto recon_man = &reconstruction_managers[cluster];
  //   std::cout << "cluster = " << cluster << ", cluster->image_ids.size() = " << cluster->image_ids.size() << std::endl;
  //   // ReconstructionManager* reclu1 = &reconstruction_managers[cluster];

    
  //   //FIXME 开始重建工作
  //   thread_pool.AddTask(ReconstructCluster, *cluster,
  //                       &reconstruction_managers[cluster]); 
  //                       //[]访问或插入元素，如果该key不存在则执行插入；反之，覆盖掉原来的key
  // }
  // thread_pool.Wait();

  //////////////////////////////////////////////////////////////////////////////
  // Merge clusters
  //////////////////////////////////////////////////////////////////////////////


  PrintHeading1("Merging clusters");

  // const colmap::SceneClustering::Cluster* scene_clusters = new colmap::SceneClustering::Cluster();

  // scene_clusters->child_clusters.push_back(*c1);
  // scene_clusters->child_clusters.push_back(*c0);
  // scene_clusters->child_clusters.push_back(*c2);
  // scene_clusters->child_clusters.push_back(*c3);
  // scene_clusters->child_clusters.push_back(*c4);

  // scene_clustering.GetRootCluster() = &scene_clusters;

  // std::cout << " scene_clusters->child_clusters.size = " << scene_clusters->child_clusters.size() << std::endl;
  // std::cout << " scene_clusters 的地址为 ：" << scene_clusters << std::endl;

  // std::cout << " scene_cluster 里面元素的地址为 ：" << std::endl;
  for(const auto& cluster : scene_clustering.GetLeafClusters()){
    std::cout << "cluster.image.size() = " << cluster.image_ids.size() << ", 地址" << &cluster << std::endl;
    const SceneClustering::Cluster * clus = new colmap::SceneClustering::Cluster();
    clus = &cluster;
    ReconstructionManager rec = reconstruction_managers[cluster];
    thread_pool.AddTask(ReconstructCluster, cluster,
                        &reconstruction_managers[&cluster]); 
  }
  thread_pool.Wait();

  std::cout << " *scene_clustering.GetRootClusters() = " << scene_clustering.GetRootCluster() << std::endl;
  std::cout << " *scene_clustering.GetRootCluster()->image_ids.size() = " << scene_clustering.GetRootCluster()->image_ids.size() << std::endl;

  std::cout << "&reconstruction_managers = " << &reconstruction_managers << std::endl;
  for (const auto& recon_man : reconstruction_managers)
  {
    std::cout << "reconstruction_managers.first(Cluster *) = " << recon_man.first << " , recon_man.first->image_ids.size() =  "  << 
    recon_man.first->image_ids.size() << " , recon_man.second.Size() = " << recon_man.second.Size() << std::endl;
  }

  // MergeClusters(*scene_clusters, &reconstruction_managers);

  MergeClusters(*scene_clustering.GetRootCluster(), &reconstruction_managers);

  // FIXME reconstruction_manager_ 为 class HierarchicalMapperController的private成员变量
  // class HierarchicalMapperController继承了thread因此HierMapper线程跑起来的时候这个自动被执行
  // std::unordered_map : begin()指向unordered_map 容器的 container(1)或者its buckets(2) 
  // std::unordered_map是无序的，因此begin()并不指向某个具体的数
  // std::unordered_map second指向mapped valued值
  CHECK_EQ(reconstruction_managers.size(), 1);
  *reconstruction_manager_ = std::move(reconstruction_managers.begin()->second);

  std::cout << std::endl;
  GetTimer().PrintMinutes();
}

}  // namespace colmap
