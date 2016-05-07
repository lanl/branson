/*
  Author: Alex Long
  Date: 5/2/201
  Name: load_balance.h
*/

#ifndef load_balance_h_
#define load_balance_h_

#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <vector>

std::unordered_map<int32_t, int32_t> load_balance(const int& rank, 
  const int& n_rank, 
  const uint64_t n_particle) 
{
  using std::unordered_map;
  using std::unordered_set;
  using std::vector;

  //////////////////////////////////////////////////////////////////////////////
  // Calculate load imbalance for each rank
  //////////////////////////////////////////////////////////////////////////////

  vector<int64_t> domain_particles(n_rank, 0);
  vector<int64_t> domain_delta_p(n_rank, 0);

  domain_particles[rank] = n_particle;
  
  MPI_Allreduce(MPI_IN_PLACE, &domain_particles[0], n_rank, MPI_UNSIGNED_LONG,
    MPI_SUM, MPI_COMM_WORLD);

  int64_t n_global_particles = 
    std::accumulate(domain_particles.begin(), domain_particles.end(),0);
  
  int64_t n_balance_particles = n_global_particles/n_rank;

  for (uint32_t ir=0; ir<n_rank; ir++) 
    domain_delta_p[ir] = domain_particles[ir] - n_balance_particles;

  //////////////////////////////////////////////////////////////////////////////
  // Determine which ranks to send or reiceve work from  
  //////////////////////////////////////////////////////////////////////////////

  unordered_map<int32_t, int32_t> n_send_rank;
  vector<uint32_t> work_donor_ranks;

  int send_to_nbr, rank_delta_p, left_node, right_node;
  for (uint32_t ir=0; ir<n_rank; ir++) {
    rank_delta_p = domain_delta_p[ir];
    left_node = ir-1;
    right_node = ir+1;
    while (rank_delta_p > 0 && (left_node>=0 || right_node<n_rank)) {
      if (left_node >= 0) {
        if (domain_delta_p[left_node] < 0 && rank_delta_p > 0) {
          send_to_nbr = abs(domain_delta_p[left_node]);
          if (send_to_nbr > rank_delta_p) send_to_nbr = rank_delta_p;

          // if this is your rank, remember it in your map
          if (rank==ir) n_send_rank[left_node] = send_to_nbr;

          // if you are going to receive work, remember donor rank
          if (rank==left_node) work_donor_ranks.push_back(ir);

          domain_delta_p[left_node]+=send_to_nbr;
          rank_delta_p -= send_to_nbr;
        }
        left_node--;
      } // if left_node >= 0
      if (right_node < n_rank) {
        if (domain_delta_p[right_node] < 0 && rank_delta_p > 0) {
          send_to_nbr = abs(domain_delta_p[right_node]);
          if (send_to_nbr > rank_delta_p) send_to_nbr = rank_delta_p;

          // if this is your rank, remember it in your map
          if (rank==ir) n_send_rank[left_node] = send_to_nbr;
          // if you are going to receive work, remember donor rank

          if (rank==left_node) work_donor_ranks.push_back(ir);

          domain_delta_p[right_node]+=send_to_nbr;
          rank_delta_p -= send_to_nbr;
        }        
        right_node++;
      } // if ir+dist < n_rank
    } //end while
  } // end for ir in n_rank

  //////////////////////////////////////////////////////////////////////////////
  // Send and receive work to balance the load
  //////////////////////////////////////////////////////////////////////////////

  // make MPI requests for the number of work sources
  uint32_t n_donors = work_donor_rank.size();
  MPI_Request *work_recv_request = new MPI_Request[n_donors];
  MPI_Request *phtn_recv_request = new MPI_Request[n_donors];
  vector<Buffer<Photon> > photon_buffer(n_donors);
  vector<Buffer<Work_Packet> > work_buffer(n_donors);
  // get the total number of particles deficient on this rank to size buffers
  uint32_t n_deficient = n_balance - n_particle;

  // post receives from ranks that are sending you work
  for (uint32_t i=0; i<n_donors; i++) {
    // make buffer maximum possible size
    work_buffer[i].resize(n_deficient);

    // post work packet receives
    MPI_Irecv(work_buffer[i].get_buffer(), n_deficient, MPI_WPacket, 
      work_donor_ranks[i], work_tag, MPI_COMM_WORLD, &work_recv_request[i]);

    // make buffer maximum possible size
    photon_buffer[i].resize(n_deficient);

    // post particle receives
    MPI_Irecv(phtn_buffer[i].get_buffer(), n_deficient, MPI_Particle, 
      work_donor_ranks[i], MPI_Particle, photon_tag, MPI_COMM_WORLD, 
      &phtn_recv_request[i]);
  }

  uint32_t n_acceptors =  n_send_rank.size();
  MPI_Request *work_send_request = new MPI_Request[n_acceptors];
  MPI_Request *phtn_send_request = new MPI_Request[n_acceptors];
  uint32_t dest_rank;
  uint32_t temp_n_send;
  uint32_t ireq;  //! Request index

  // iterate over unordered map and prepare work for other ranks
  for (std::unordered_map<int32_t,int32_t>::iterator map_itr=n_send_rank.begin();
    map_itr !=n_send_rank.end(); map_itr++) {

    dest_rank = map_itr->first;
    temp_n_send = map_itr->second;

    vector<Work_Packet> work_to_send;
    Work_Packet temp_packet, leftover_packet;
    // first, try to send work packets to the other ranks
    while (!work_packets_stack.empty() && temp_n_send > 0) {
      temp_packet = work_packet_stack.top();
      // split work if it's larger than the number needed
      if (temp_packet.get_n_particles() > temp_n_send) {
        leftover_packet = work_packet.split(temp_n_send);
        // pop off the original packet
        work_packet_stack.pop();
        // push on the leftover packet
        work_packet_stack.push(leftover_packet); 
      }
      // otherwise, pop the temp work packet off the stack 
      else work_packet_stack.pop();

      // add packet to send list
      work_to_send.push_back(temp_packet);
      // subtract particle in packet from temp_n_send
      temp_n_send -= temp_packet.get_n_particles();
    }

    // send census particles instead (not preferred, these photons are likely
    // to travel farther and thus require more memory
    vector<Photon> photons_to_send;
    while (!census_list.empty() && temp_n_send > 0) {
      if (census_list.size() > temp_n_send) {
        photons_to_send.insert(photons_to_send.begin(),
          census_list.begin(), census_list.begin()+temp_n_send);
      }
    }

    // send both work and photon vectors, even if they're empty
    MPI_Isend(&work_to_send[0], work_to_send.size(), MPI_WPacket,
      dest_rank, work_tag, MPI_COMM_WORLD, &work_send_request[ireq]);
    MPI_Isend(&photons_to_send[0], photons_to_send.size(), MPI_Particle,
      dest_rank, photon_tag, MPI_COMM_WORLD, &phtn_send_request[ireq]);
    // increment the request index counter
    ireq++; 
  }
}  
#endif // load_balance_h_
