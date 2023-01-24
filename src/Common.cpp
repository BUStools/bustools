#include "Common.hpp"
#include <fstream>

std::vector<int32_t> intersect(std::vector<int32_t> &u, std::vector<int32_t> &v) {
  std::vector<int32_t> res;
  auto a = u.begin();
  auto b = v.begin();

  while (a != u.end() && b != v.end()) {
    if (*a < *b) {
      ++a;
    } else if (*b < *a) {
      ++b;
    } else {
      // match
      res.push_back(*a);
      ++a;
      ++b;
    }
  }
  return std::move(res);
}

std::vector<int32_t> union_vectors(const std::vector<std::vector<int32_t>> &v) {
  std::vector<int32_t> u;
  for (const auto &vv : v) {
    for (auto t : vv) {
      u.push_back(t);
    }
  }
  
  std::sort(u.begin(), u.end());
  u.erase(std::unique(u.begin(), u.end()), u.end());
  return std::move(u);
}

std::vector<int32_t> intersect_vectors(const std::vector<std::vector<int32_t>> &v) {
  std::vector<int32_t> u;

  if (!v.empty()) {
    u = v[0]; // copy
    for (size_t i = 1; i < v.size(); i++) {
      auto &vv = v[i];
      int j = 0;
      int k = 0;
      int l = 0;
      int n = u.size();
      int m = vv.size();
      // u and v are sorted, j,k,l = 0
      while (j < n && l < m) {
        // invariant: u[:k] is the intersection of u[:j] and vv[:l], j <= n, l <= m
        //            u[:j] <= u[j:], vv[:l] <= v[l:], u[j:] is sorted, v[l:] is sorted, u[:k] is sorted
        if (u[j] < vv[l]) {
          j++;
        } else if (u[j] > vv[l]) {
          l++;
        } else {
          // match
          if (k < j) {
            std::swap(u[k], u[j]);
          }
          k++;
          j++;
          l++;
        }
      }
      if (k < n) {
        u.resize(k);
      }    
    }
  }

  return std::move(u);
}

int32_t intersect_ecs(const std::vector<int32_t> &ecs, std::vector<int32_t> &u, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ecmap, u_map_<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv, std::vector<std::vector<int32_t>> &ec2genes) {
  if (ecs.empty()) {
    return -1;
  }

  if (ecs[0] < 0 || ecs[0] >= ecmap.size()) {
    return -1;
  }
  
  if (ecs.size() == 1) {
    return ecs[0]; // no work
  }

  u.resize(0);
  auto &v = ecmap[ecs[0]]; // copy
  for (size_t i = 0; i< v.size(); i++) {
    u.push_back(v[i]);
  }

  for (size_t i = 1; i < ecs.size(); i++) {
    if (ecs[i] < 0 || ecs[i] >= ecmap.size()) {
      return -1;
    }
    const auto &v = ecmap[ecs[i]];

    int j = 0;
    int k = 0;
    int l = 0;
    int n = u.size();
    int m = v.size();
    // u and v are sorted, j,k,l = 0
    while (j < n && l < m) {
      // invariant: u[:k] is the intersection of u[:j] and v[:l], j <= n, l <= m
      //            u[:j] <= u[j:], v[:l] <= v[l:], u[j:] is sorted, v[l:] is sorted, u[:k] is sorted
      if (u[j] < v[l]) {
        j++;
      } else if (u[j] > v[l]) {
        l++;
      } else {
        // match
        if (k < j) {
          std::swap(u[k], u[j]);
        }
        k++;
        j++;
        l++;
      }
    }
    if (k < n) {
      u.resize(k);
    }
  }

  if (u.empty()) {
    return -1;
  }
  auto iit = ecmapinv.find(u);
  if (iit == ecmapinv.end()) { 
    // create new equivalence class
    int32_t ec = ecmap.size();
    ecmap.push_back(u);
    ecmapinv.insert({u,ec});
    // figure out the gene list
    std::vector<int32_t> v;
    vt2gene(u, genemap, v);
    ec2genes.push_back(std::move(v));
    return ec;
  } else {
    return iit->second;
  }
}


void vt2gene(const std::vector<int32_t> &v, const std::vector<int32_t> &genemap, std::vector<int32_t> &glist) {
  int lastg = -2;
  int n = v.size();

  for (int i = 0; i < n; i++) {
    auto t = v[i];
    auto g = genemap[t];

    if (g != lastg && g != -1) {
      glist.push_back(g);
      lastg = g;
    }
  }

  if (glist.size() > 1) {
    // sort and remove duplicates
    std::sort(glist.begin(), glist.end());
    glist.erase(std::unique(glist.begin(), glist.end()), glist.end());
  }
}


void intersect_genes_of_ecs(const std::vector<int32_t> &ecs, const  std::vector<std::vector<int32_t>> &ec2genes, std::vector<int32_t> &glist) {
  glist.resize(0);
  if  (ecs.empty()) {
    return;
  }
  // copy first to glist
  const auto &v = ec2genes[ecs[0]];
  for (auto x : v) {
    glist.push_back(x);
  }
  // intersect the rest
  for (int i = 1; i < ecs.size(); i++) {
    const auto &v = ec2genes[ecs[i]];

    int j = 0;
    int k = 0;
    int l = 0;
    int n = glist.size();
    int m = v.size();
    // u and v are sorted, j,k,l = 0
    while (j < n && l < m) {
      // invariant: u[:k] is the intersection of u[:j] and v[:l], j <= n, l <= m
      //            u[:j] <= u[j:], v[:l] <= v[l:], u[j:] is sorted, v[l:] is sorted, u[:k] is sorted
      if (glist[j] < v[l]) {
        j++;
      } else if (glist[j] > v[l]) {
        l++;
      } else {
        // match
        if (k < j) {
          std::swap(glist[k], glist[j]);
        }
        k++;
        j++;
        l++;
      }
    }
    if (k < n) {
      glist.resize(k);
    }
  }
}


int32_t intersect_ecs_with_genes(const std::vector<int32_t> &ecs, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ecmap, u_map_<std::vector<int32_t>, int32_t, SortedVectorHasher> &ecmapinv, std::vector<std::vector<int32_t>> &ec2genes, bool assumeIntersectionIsEmpty) {
  
  std::vector<std::vector<int32_t>> gu; // per gene transcript results
  std::vector<int32_t> u; // final list of transcripts
  std::vector<int32_t> glist;

  int32_t lastg = -2;
  // todo, replace by intersection of the genelist
  for (const auto ec : ecs) {
    auto g = ec2genes[ec];
    if (g.size() == 1 && g[0] != lastg) {
      glist.push_back(g[0]);
      lastg = g[0];
    } else if (g.size() > 1) {
      lastg = -2;
      for (auto &x : g) {
        glist.push_back(x);
      }
    }
  }
  
  if (glist.empty()) {
    return -1;
  }

  // sort and remove unique
  std::sort(glist.begin(), glist.end());
  glist.erase(std::unique(glist.begin(), glist.end()), glist.end());

  if (glist.size() == 1 && assumeIntersectionIsEmpty) {
    // frequent case, single gene replace with union
    for (auto ec : ecs) {
      for (const auto &t : ecmap[ec]) {      
        u.push_back(t);
      }
    }
    std::sort(u.begin(), u.end());
    u.erase(std::unique(u.begin(), u.end()), u.end());

    // look up ecs based on u
    int32_t ec = -1;
    
    auto it = ecmapinv.find(u);
    if (it != ecmapinv.end()) {
      ec = it->second;              
    } else {
      ec = ecmapinv.size();
      ecmapinv.insert({u,ec});  
      ecmap.push_back(u);
      std::vector<int32_t> v;
      vt2gene(u, genemap, v);
      ec2genes.push_back(std::move(v));
    }

    return ec; // done
  } else {
    // separate per gene
    for (auto g : glist) {
      gu.clear();
      
      for (auto ec : ecs) {
        std::vector<int32_t> tg;
        for (const auto &t : ecmap[ec]) {
          if (genemap[t] == g) {
            tg.push_back(t);
          }
        }
        if (!tg.empty()) {
          gu.push_back(std::move(tg));
        }
      }
      auto uu = intersect_vectors(gu);

      // if gene intersection is empty, use union
      if (uu.empty()) {
        uu = union_vectors(gu);
      }

      for (auto t : uu) { 
        u.push_back(t);
      }
    }

    if (u.empty()) {
      return -1;
    }
    std::sort(u.begin(), u.end());

    int32_t ec = -1;
    auto it = ecmapinv.find(u);
    if (it != ecmapinv.end()) {
      ec = it->second;              
    } else {
      ec = ecmapinv.size();
      ecmapinv.insert({u,ec});
      ecmap.push_back(u);
      std::vector<int32_t> v;
      vt2gene(u, genemap, v);
      ec2genes.push_back(std::move(v));
    }
    return ec;
  } 
  
}


void create_ec2genes(const std::vector<std::vector<int32_t>> &ecmap, const std::vector<int32_t> &genemap, std::vector<std::vector<int32_t>> &ec2gene) {
  std::vector<int32_t> u;
  for (int ec = 0; ec < ecmap.size(); ec++) {
    const auto &v = ecmap[ec];    
    vt2gene(v, genemap, u);
    ec2gene.push_back(std::move(u));
  }
}

COUNT_MTX_TYPE intersect_ecs_with_subset_txs(int32_t ec, const std::vector<std::vector<int32_t>> &ecmap, const std::vector<int32_t>& tx_split) {
  if (tx_split.size() == 0) return COUNT_DEFAULT;
  std::vector<int32_t> ecs;
  ecs.push_back(ec);
  return intersect_ecs_with_subset_txs(ecs, ecmap, tx_split);
}

COUNT_MTX_TYPE intersect_ecs_with_subset_txs(const std::vector<int32_t>& ecs, const std::vector<std::vector<int32_t>> &ecmap, const std::vector<int32_t>& tx_split) {
  if (tx_split.size() == 0) return COUNT_DEFAULT;
  if (ecs.size() == 0) return COUNT_AMBIGUOUS; // Shouldn't happen
  size_t n_1 = 0;
  size_t n_2 = 0;
  for (auto ec : ecs) { // We still need to optimize this
    for (auto t: ecmap[ec]) {
      if(std::find(tx_split.begin(), tx_split.end(), t) != tx_split.end()) {
        n_2++;
      } else {
        n_1++;
      }
      if (n_1 > 0 && n_2 > 0) break; // Stop searching
    }
    if (n_1 > 0 && n_2 > 0) break; // Stop searching
  }
  return (n_1 > 0 && n_2 > 0 ? COUNT_AMBIGUOUS : (n_1 > 0 ? COUNT_DEFAULT : COUNT_SPLIT));
}


void copy_file(std::string src, std::string dest) {
	std::ifstream  isrc(src, std::ios::binary);
	std::ofstream  idest(dest, std::ios::binary);

	idest << isrc.rdbuf();
}