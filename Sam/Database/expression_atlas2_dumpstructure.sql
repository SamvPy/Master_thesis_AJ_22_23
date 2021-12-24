-- MySQL dump 10.13  Distrib 8.0.22, for Win64 (x86_64)
--
-- Host: localhost    Database: expression_atlas2
-- ------------------------------------------------------
-- Server version	8.0.22

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `assay`
--
USE expression_atlas_shotgun;
DROP TABLE IF EXISTS `assay`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `assay` (
  `assay_id` bigint NOT NULL AUTO_INCREMENT,
  `project_id` int NOT NULL,
  `filename` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`assay_id`,`project_id`),
  UNIQUE KEY `id_UNIQUE` (`assay_id`),
  KEY `fk_assay_project1_idx` (`project_id`),
  CONSTRAINT `fk_assay_project1` FOREIGN KEY (`project_id`) REFERENCES `project` (`project_id`)
) ENGINE=InnoDB AUTO_INCREMENT=27704 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `assay_modifications`
--

DROP TABLE IF EXISTS `assay_modifications`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `assay_modifications` (
  `assay_id` bigint NOT NULL,
  `mod_id` int NOT NULL,
  `count` int DEFAULT '1',
  UNIQUE KEY `uq_assay_modifications` (`assay_id`,`mod_id`),
  KEY `assay_modifications_ibfk_2` (`mod_id`),
  CONSTRAINT `assay_modifications_ibfk_1` FOREIGN KEY (`assay_id`) REFERENCES `assay` (`assay_id`),
  CONSTRAINT `assay_modifications_ibfk_2` FOREIGN KEY (`mod_id`) REFERENCES `modifications` (`mod_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene_ontology`
--

DROP TABLE IF EXISTS `gene_ontology`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `gene_ontology` (
  `go_id` varchar(20) NOT NULL,
  `go_term` varchar(255) DEFAULT NULL,
  `go_domain` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`go_id`),
  UNIQUE KEY `id_UNIQUE` (`go_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mod_freq`
--

DROP TABLE IF EXISTS `mod_freq`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `mod_freq` (
  `mod_id` int NOT NULL,
  `freq` bigint NOT NULL DEFAULT '0',
  `mod_type` varchar(255) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `modifications`
--

DROP TABLE IF EXISTS `modifications`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `modifications` (
  `mod_id` int NOT NULL AUTO_INCREMENT,
  `modification_type` varchar(255) DEFAULT NULL,
  `mass_difference` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`mod_id`),
  UNIQUE KEY `id_UNIQUE` (`mod_id`)
) ENGINE=InnoDB AUTO_INCREMENT=2030 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `peptide`
--

DROP TABLE IF EXISTS `peptide`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `peptide` (
  `peptide_id` bigint NOT NULL AUTO_INCREMENT,
  `peptide_sequence` varchar(255) DEFAULT NULL,
  `start_position` int DEFAULT NULL,
  PRIMARY KEY (`peptide_id`),
  UNIQUE KEY `id_UNIQUE` (`peptide_id`),
  UNIQUE KEY `peptide_sequence_UNIQUE` (`peptide_sequence`)
) ENGINE=InnoDB AUTO_INCREMENT=104508336 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `peptide_modifications`
--

DROP TABLE IF EXISTS `peptide_modifications`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `peptide_modifications` (
  `peptide_id` bigint NOT NULL,
  `location` varchar(20) DEFAULT NULL,
  `mod_id` int NOT NULL,
  `assay_id` bigint NOT NULL,
  UNIQUE KEY `uq_peptide_modifications` (`peptide_id`,`location`,`mod_id`,`assay_id`),
  KEY `peptide_modifications_ibfk_2` (`mod_id`),
  KEY `peptide_modifications_ibfk_3` (`assay_id`),
  CONSTRAINT `peptide_modifications_ibfk_1` FOREIGN KEY (`peptide_id`) REFERENCES `peptide` (`peptide_id`),
  CONSTRAINT `peptide_modifications_ibfk_2` FOREIGN KEY (`mod_id`) REFERENCES `modifications` (`mod_id`),
  CONSTRAINT `peptide_modifications_ibfk_3` FOREIGN KEY (`assay_id`) REFERENCES `assay` (`assay_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `peptide_to_assay`
--

DROP TABLE IF EXISTS `peptide_to_assay`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `peptide_to_assay` (
  `peptide_id` bigint NOT NULL,
  `assay_id` bigint NOT NULL,
  `quantification` float NOT NULL,
  UNIQUE KEY `uq_peptide_to_assay` (`peptide_id`,`assay_id`),
  KEY `peptide_to_assay_ibfk_2` (`assay_id`),
  CONSTRAINT `peptide_to_assay_ibfk_1` FOREIGN KEY (`peptide_id`) REFERENCES `peptide` (`peptide_id`),
  CONSTRAINT `peptide_to_assay_ibfk_2` FOREIGN KEY (`assay_id`) REFERENCES `assay` (`assay_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `peptide_to_protein`
--

DROP TABLE IF EXISTS `peptide_to_protein`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `peptide_to_protein` (
  `peptide_id` bigint NOT NULL,
  `uniprot_id` varchar(200) NOT NULL,
  UNIQUE KEY `uq_peptide_to_protein` (`peptide_id`,`uniprot_id`),
  KEY `peptide_to_protein_ibfk_2` (`uniprot_id`),
  CONSTRAINT `peptide_to_protein_ibfk_1` FOREIGN KEY (`peptide_id`) REFERENCES `peptide` (`peptide_id`),
  CONSTRAINT `peptide_to_protein_ibfk_2` FOREIGN KEY (`uniprot_id`) REFERENCES `protein` (`uniprot_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `project`
--

DROP TABLE IF EXISTS `project`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `project` (
  `project_id` int NOT NULL AUTO_INCREMENT,
  `PXD_accession` varchar(45) DEFAULT NULL,
  `experiment_type` varchar(255) DEFAULT NULL,
  `instrument` varchar(255) DEFAULT NULL,
  `keywords` longtext CHARACTER SET utf8 COLLATE utf8_general_ci NOT NULL,
  `ref` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`project_id`, `PXD_accession`)
) ENGINE=InnoDB AUTO_INCREMENT=1701 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein`
--

DROP TABLE IF EXISTS `protein`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `protein` (
  `uniprot_id` varchar(255) NOT NULL,
  `description` longtext,
  `length` varchar(255) DEFAULT NULL,
  `sequence` longtext,
  PRIMARY KEY (`uniprot_id`),
  UNIQUE KEY `id_UNIQUE` (`uniprot_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_protein_interactions`
--

DROP TABLE IF EXISTS `protein_protein_interactions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `protein_protein_interactions` (
  `ppi_id` varchar(20) NOT NULL,
  `protein1_id` varchar(20) NOT NULL,
  `protein2_id` varchar(20) NOT NULL,
  PRIMARY KEY (`ppi_id`),
  UNIQUE KEY `interaction_id_UNIQUE` (`ppi_id`),
  KEY `protein_protein_interactions_ibfk_2_idx` (`protein2_id`),
  KEY `protein_protein_interactions_ibfk_1_idx` (`protein1_id`),
  CONSTRAINT `protein_protein_interactions_ibfk_1` FOREIGN KEY (`protein1_id`) REFERENCES `protein` (`uniprot_id`),
  CONSTRAINT `protein_protein_interactions_ibfk_2` FOREIGN KEY (`protein2_id`) REFERENCES `protein` (`uniprot_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_to_go`
--

DROP TABLE IF EXISTS `protein_to_go`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `protein_to_go` (
  `uniprot_id` varchar(200) NOT NULL,
  `go_id` varchar(20) NOT NULL,
  UNIQUE KEY `uq_protein_to_go` (`uniprot_id`,`go_id`),
  KEY `protein_to_go_ibfk_2` (`go_id`),
  CONSTRAINT `protein_to_go_ibfk_1` FOREIGN KEY (`uniprot_id`) REFERENCES `protein` (`uniprot_id`),
  CONSTRAINT `protein_to_go_ibfk_2` FOREIGN KEY (`go_id`) REFERENCES `gene_ontology` (`go_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `psm`
--

DROP TABLE IF EXISTS `psm`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `psm` (
  `spectrum_id` bigint NOT NULL,
  `peptide_id` bigint NOT NULL,
  KEY `psm_ibfk_1` (`spectrum_id`),
  KEY `psm_ibfk_2` (`peptide_id`),
  CONSTRAINT `psm_ibfk_1` FOREIGN KEY (`spectrum_id`) REFERENCES `spectra` (`spectrum_id`),
  CONSTRAINT `psm_ibfk_2` FOREIGN KEY (`peptide_id`) REFERENCES `peptide` (`peptide_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `spectra`
--

DROP TABLE IF EXISTS `spectra`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `spectra` (
  `spectrum_id` bigint NOT NULL AUTO_INCREMENT,
  `mass_to_charge_ratio` int DEFAULT NULL,
  `charge` int DEFAULT NULL,
  `intensity` int DEFAULT NULL,
  PRIMARY KEY (`spectrum_id`),
  UNIQUE KEY `id_UNIQUE` (`spectrum_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tissue`
--

DROP TABLE IF EXISTS `tissue`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `tissue` (
  `tissue_id` int NOT NULL AUTO_INCREMENT,
  `tissue_name` varchar(255) DEFAULT NULL,
  `cell_type` varchar(255) DEFAULT NULL,
  `disease_status` varchar(255) DEFAULT NULL,
  `organ_id` int DEFAULT NULL,
  PRIMARY KEY (`tissue_id`),
  KEY `organ_id` (`organ_id`)
) ENGINE=InnoDB AUTO_INCREMENT=452 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tissue_to_assay`
--

DROP TABLE IF EXISTS `tissue_to_assay`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `tissue_to_assay` (
  `assay_id` bigint NOT NULL,
  `tissue_id` int NOT NULL,
  UNIQUE KEY `uq_tissue_to_assay` (`assay_id`,`tissue_id`),
  KEY `tissue_to_assay_ibfk_2` (`tissue_id`),
  CONSTRAINT `tissue_to_assay_ibfk_1` FOREIGN KEY (`assay_id`) REFERENCES `assay` (`assay_id`),
  CONSTRAINT `tissue_to_assay_ibfk_2` FOREIGN KEY (`tissue_id`) REFERENCES `tissue` (`tissue_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tissue_to_protein`
--

DROP TABLE IF EXISTS `tissue_to_protein`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `tissue_to_protein` (
  `uniprot_id` varchar(200) NOT NULL,
  `organ_id` int NOT NULL,
  `tissue` varchar(128) NOT NULL,
  `cell_type` varchar(128) NOT NULL,
  `level` int NOT NULL,
  UNIQUE KEY `uq_tissue_to_protein` (`uniprot_id`,`organ_id`),
  KEY `tissue_to_protein_ibfk_2` (`organ_id`),
  CONSTRAINT `tissue_to_protein_ibfk_1` FOREIGN KEY (`uniprot_id`) REFERENCES `protein` (`uniprot_id`),
  CONSTRAINT `tissue_to_protein_ibfk_2` FOREIGN KEY (`organ_id`) REFERENCES `tissue` (`tissue_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2021-09-23 13:47:22
