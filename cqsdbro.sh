sudo umount /mnt/disks/cqsdb
gcloud compute instances detach-disk cqscluster-controller --disk=cqsdb
gcloud compute instances attach-disk cqscluster-controller --disk=cqsdb --mode=ro
sudo mount -o discard,defaults /dev/sdb /mnt/disks/cqsdb
