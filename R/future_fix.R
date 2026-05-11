# 1. 提高future包的全局对象大小限制
# 添加在脚本开始或使用future前

# 导入必要的包
library(future)
library(future.apply)

# 方法1: 增加全局对象大小限制（临时解决方案）
options(future.globals.maxSize = 8000 * 1024^2)  # 设置为8GB

# 方法2: 优化future函数调用（推荐长期解决方案）
# 创建一个优化版本的future_lapply函数
optimized_future_lapply <- function(X, FUN, ..., future.seed = TRUE) {
  # 1. 减少传递给worker的对象大小
  # 如果FUN是大型函数或引用了大型数据结构
  FUN_optimized <- function(...) {
    # 在worker中重新创建需要的数据而不是传递
    # 示例: 如果需要在worker中使用某个大型对象
    # large_data <- readRDS("path/to/large_data.rds") 
    
    # 调用原始函数
    FUN(...)
  }
  
  # 2. 使用future.apply::future_lapply进行并行处理
  future.apply::future_lapply(X, FUN_optimized, ..., future.seed = future.seed)
}

# 示例改进代码 - 用于处理大型Seurat对象的并行函数
process_large_seurat_parallel <- function(seurat_obj, chunk_ids, process_fn) {
  # 确保worker不需要整个seurat对象
  # 只提取必要的数据给worker
  
  # 按照cell chunk创建处理任务
  cell_chunks <- split(colnames(seurat_obj), cut(seq_along(colnames(seurat_obj)), chunk_ids))
  
  # 使用优化版本并行处理，每个worker处理一组细胞
  results <- optimized_future_lapply(cell_chunks, function(cells) {
    # 这里不要传递整个seurat_obj
    # 只提取必要的数据
    required_data <- list(
      counts = as.matrix(seurat_obj[["RNA"]]@counts[, cells]),
      metadata = seurat_obj@meta.data[cells, ]
    )
    
    # 处理提取的数据
    process_fn(required_data)
  })
  
  # 合并结果
  combined_results <- do.call(rbind, results)
  return(combined_results)
}

# 使用方法示例
run_parallel_analysis <- function(seurat_obj) {
  # 设置并行后端
  plan(multisession, workers = 4)  # 根据您的系统调整
  
  # 将分析分成多个chunks
  n_chunks <- 10
  
  # 使用优化的并行处理函数
  results <- process_large_seurat_parallel(
    seurat_obj = seurat_obj,
    chunk_ids = n_chunks,
    process_fn = function(data) {
      # 在这里处理每个chunk的数据
      # 不要在这里引用大型全局对象
      return(data$metadata)
    }
  )
  
  # 恢复顺序执行
  plan(sequential)
  
  return(results)
}
