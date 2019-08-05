#pragma once
class PCLViewer
{
public:
	PCLViewer();
	~PCLViewer();

	void ReadPcdFile(pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud, const fs::path &fileName);		//	��ȡpcd�ļ��еõ�
	void AscToPcd(const fs::path &fileName);														//	��asc�ļ�ת��pcd�ļ���ʽ
private:
	void ReadBuffer(const fs::path &fileDir, char **buffer);
private:
	unsigned int red;
	unsigned int green;
	unsigned int blue;
};

